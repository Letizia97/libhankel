"""
Speed of the SESANS Hankel transform: sasmodels' matrix against LibHankel.

The two do the same integral in very different ways, and the costs are not
directly comparable without saying what is being timed:

* sasmodels pays a large one-off cost to build the SesansTransform matrix H
  (its size is set by the spin-echo lengths, not by the model), then each model
  evaluation is I(q) on the q grid plus a matrix-vector product. In a fit the
  build happens once and the evaluation happens once per iteration.
* LibHankel has no setup at all; the whole cost is in the call, and it scales
  with the number of spin-echo lengths and with the requested tolerance.

So the fair comparison is LibHankel's call against sasmodels' *per evaluation*
cost, with the build cost reported separately. Accuracy is reported alongside,
against the analytic SESANS correlation function of a sphere, because a
transform is only as useful as the tolerance it holds.

Run with:

    python examples/sesans/benchmark.py
"""

import time

import libhankel
import numpy as np
from sasmodels.core import load_model
from sasmodels.data import empty_sesans
from sasmodels.direct_model import DirectModel

RADIUS = 1000.0  # A
CONTRAST = 1.0
WAVELENGTH = 5.0  # A
N_XI = 50
# Best-of-N, not an average: on a laptop the slow runs are scheduling noise,
# and N has to be biggish for the sasmodels evaluation to settle down.
REPEATS = 15

MODEL_PARS = dict(radius=RADIUS, sld=CONTRAST, sld_solvent=0.0, scale=1.0, background=0.0)

# G(0) for a sphere, with the 1 / 2 pi convention used by SesansTransform.
G_ZERO = 2 * np.pi * CONTRAST**2 * RADIUS**4

# n_eval / eps_rel are set to values that converge on both grids below; the
# QWE strategies fail on the linear grid at the default n_eval of 250, because
# G(xi) passes through zero at xi = 2 R and a relative tolerance cannot be met
# there. f_max for Ogata has to be small for this integrand.
STRATEGIES = [
    ("DHT_Guptasarma_Fast", {}),
    ("DHT_Guptasarma", {}),
    ("DHT_Key_51", {}),
    ("DHT_Key_101", {}),
    ("DHT_Key_201", {}),
    ("DHT_Anderson_801", {}),
    ("Fixed_DE_Ogata", {"n_eval": 250, "f_max": 1e-4}),
    ("Adaptive_DE_Ooura", {"n_eval": 250, "eps_rel": 1e-9}),
    ("QWE_Key", {"n_eval": 500, "eps_rel": 1e-8}),
    ("QWE_Chave", {"n_eval": 500, "eps_rel": 1e-8}),
]


def analytic_sesans_sphere(xi):
    """
    P(xi) / G(0) for a homogeneous sphere, the reference both routes approximate.

    Same expression as ``compute_analytical_spheres`` in
    ``src/utils/analytical_form_factors.c``, normalised by G(0) = 2 pi eta^2 R^4.
    """
    u = np.clip(xi / (2.0 * RADIUS), 0.0, 1.0)
    root = np.sqrt(1.0 - u**2)
    log_term = np.log(np.maximum(u / (1.0 + root), 1e-15))
    return (root * (2.0 + u**2) + u**2 * (4.0 - u**2) * log_term - 2.0) / 2.0


def timings(call, repeats=REPEATS):
    """Median and best wall-clock time over *repeats* calls, plus the result."""
    measured = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = call()
        measured.append(time.perf_counter() - start)
    return float(np.median(measured)), min(measured), result


def benchmark_grid(label, xi, model):
    reference = analytic_sesans_sphere(xi)

    build_time, _, _ = timings(
        lambda: DirectModel(empty_sesans(z=xi, wavelength=WAVELENGTH), model), repeats=3
    )
    calculator = DirectModel(empty_sesans(z=xi, wavelength=WAVELENGTH), model)
    calculator(**MODEL_PARS)  # warm up: first call loads and compiles the kernel
    eval_time, eval_best, polarisation = timings(lambda: calculator(**MODEL_PARS))

    n_q = len(calculator.resolution.q_calc)
    g_zero = np.dot(calculator.resolution._H0, calculator.Iq_calc)
    sasmodels_error = np.max(np.abs(polarisation / g_zero - reference))

    print(f"\n{label} spin-echo lengths, n = {len(xi)}")
    print(f"  sasmodels  build SesansTransform  {build_time * 1e3:7.2f} ms  ({n_q} q points)")
    print(
        f"  sasmodels  per model evaluation   {eval_time * 1e3:7.2f} ms"
        f"  (best {eval_best * 1e3:.2f})   max error {sasmodels_error:.1e}"
    )
    print(
        f"  {'libhankel strategy':<22} {'median':>9} {'best':>9}"
        f" {'vs sasmodels':>13} {'max error':>11}"
    )

    for name, strategy_params in STRATEGIES:
        try:
            call_time, call_best, transform = timings(
                lambda: libhankel.hankel_transform(
                    0,
                    "sphere",
                    xi,
                    np.array([RADIUS, CONTRAST]),
                    name,  # noqa: B023
                    strategy_params,  # noqa: B023
                )
            )
        except RuntimeError as error:
            print(f"  {name:<22} {'-':>9} {'-':>9} {'-':>13}   {error}")
            continue

        error = np.max(np.abs(np.asarray(transform) / (2 * np.pi) / G_ZERO - 1.0 - reference))
        print(
            f"  {name:<22} {call_time * 1e3:6.2f} ms {call_best * 1e3:6.2f} ms"
            f" {call_time / eval_time:12.2f}x {error:11.1e}"
        )


def main():
    model = load_model("sphere")
    benchmark_grid("log-spaced", np.logspace(np.log10(50.0), np.log10(5000.0), N_XI), model)
    benchmark_grid("linear", np.linspace(50.0, 5000.0, N_XI), model)


if __name__ == "__main__":
    main()
