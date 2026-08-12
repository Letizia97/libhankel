"""
The same speed comparison as ``benchmark.py``, driven by real SESANS data.

``benchmark.py`` invents its spin-echo lengths. Real measurements do not look
like that, and the shape of the data is exactly what sets the cost of the
sasmodels route: the SesansTransform picks its q grid from the spin-echo
lengths in the file (q_min from the longest, q_max from the spacing), and the
analyser acceptance and wavelength decide how much of the resulting matrix is
thrown away again.

The measured data supplies the spin-echo lengths, the wavelengths and the
acceptance angle. It does not supply I(q): the transform is always applied to
model I(q) evaluated on the grid the transform chose, and the measured
polarisation is the fit target, not an input to the transform. So this
benchmark reads real files, builds the transform they imply, and puts a sphere
model through it.

The three datasets are the ones shipped with ``sasdata``, so there is nothing
to download:

* ``sphere_isis.ses``     PMMA hard spheres in deuterated decalin, measured at
                          ISIS. Time-of-flight, so every point has its own
                          wavelength and the acceptance mask actually bites.
* ``sphere2micron.ses``   2 um polystyrene in H2O/D2O, measured at Delft.
                          Monochromatic, nothing masked.
* ``sphere2micron_long.ses``  the same sample out to 20 um.

A strategy that cannot do every spin-echo length in the file gets no result at
all, because ``hankel_transform`` raises for the whole array. Where that
happens the failing points are counted one at a time and reported, since which
xi they are is the interesting part. Those runs also write "did not converge"
to stderr from the C library, ahead of the tables.

Requires ``sasmodels`` and ``sasdata``:

    python -m pip install sasmodels sasdata
    python examples/sesans/benchmark_real_data.py
"""

import time

import libhankel
import numpy as np
from sasdata import data_path
from sasmodels.core import load_model
from sasmodels.data import load_data
from sasmodels.direct_model import DirectModel

SESANS_DIR = data_path / "sesans_data"

# Radii are the nominal particle sizes quoted for each sample. The contrast and
# the scale are nominal too: the comparison below is on P(xi) / G(0), which is
# insensitive to both, so only the radius matters for the shape of the curve.
DATASETS = [
    ("sphere_isis.ses", 1000.0, "PMMA in deuterated decalin, ISIS, time-of-flight"),
    ("sphere2micron.ses", 10000.0, "2 um polystyrene, Delft, monochromatic"),
    ("sphere2micron_long.ses", 10000.0, "the same sample, out to 20 um"),
]

CONTRAST = 1.0
MODEL_PARS = dict(sld=1.41, sld_solvent=2.70, scale=1.0, background=0.0)

# Best-of/median-of N: the sasmodels evaluation is memory-bound and scatters,
# so a single timing is not worth much.
REPEATS = 15

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


def analytic_sesans_sphere(xi, radius):
    """
    P(xi) / G(0) for a homogeneous sphere, the reference both routes approximate.

    Same expression as ``compute_analytical_spheres`` in
    ``src/utils/analytical_form_factors.c``, normalised by G(0) = 2 pi eta^2 R^4.
    """
    u = np.clip(xi / (2.0 * radius), 0.0, 1.0)
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


def count_failures(xi, radius, name, strategy_params):
    """How many individual spin-echo lengths a strategy cannot do."""
    failed = []
    for z in xi:
        try:
            libhankel.hankel_transform(
                0, "sphere", np.array([z]), np.array([radius, CONTRAST]), name, strategy_params
            )
        except RuntimeError:
            failed.append(z)
    return failed


def benchmark_dataset(filename, radius, description, model):
    data = load_data(str(SESANS_DIR / filename))
    xi = data.x
    wavelength = np.asarray(data.source.wavelength)
    theta_max, theta_unit = data.sample.zacceptance
    reference = analytic_sesans_sphere(xi, radius)

    build_time, _, _ = timings(lambda: DirectModel(data, model), repeats=3)
    calculator = DirectModel(data, model)
    calculator(radius=radius, **MODEL_PARS)  # warm up: first call compiles the kernel
    eval_time, eval_best, polarisation = timings(lambda: calculator(radius=radius, **MODEL_PARS))

    transform = calculator.resolution
    hankel_matrix = transform._H
    masked = 100.0 * (hankel_matrix == 0).sum() / hankel_matrix.size
    g_zero = np.dot(transform._H0, calculator.Iq_calc)
    sasmodels_error = np.max(np.abs(polarisation / g_zero - reference))

    print(f"\n{filename}   {description}")
    print(
        f"  {len(xi)} spin-echo lengths, {xi.min():.0f} to {xi.max():.0f} A;"
        f"  lambda {wavelength.min():.2f} to {wavelength.max():.2f} A;"
        f"  theta_max {theta_max:g} {theta_unit};  sphere radius {radius:g} A"
    )
    print(
        f"  sasmodels  build SesansTransform  {build_time * 1e3:7.2f} ms"
        f"  ({len(transform.q_calc)} q points, {hankel_matrix.nbytes / 1e6:.1f} MB,"
        f" {masked:.1f}% of H zeroed by the acceptance)"
    )
    print(
        f"  sasmodels  per model evaluation   {eval_time * 1e3:7.2f} ms"
        f"  (best {eval_best * 1e3:.2f})   max error {sasmodels_error:.1e}"
    )
    print(
        f"  {'libhankel strategy':<22} {'median':>9} {'best':>9}"
        f" {'vs sasmodels':>13} {'max error':>11}"
    )

    g_zero_analytic = 2 * np.pi * CONTRAST**2 * radius**4
    for name, strategy_params in STRATEGIES:
        try:
            call_time, call_best, result = timings(
                lambda: libhankel.hankel_transform(
                    0,
                    "sphere",
                    xi,
                    np.array([radius, CONTRAST]),
                    name,  # noqa: B023
                    strategy_params,  # noqa: B023
                )
            )
        except RuntimeError:
            failed = count_failures(xi, radius, name, strategy_params)
            near = ", ".join(f"{z:.0f}" for z in failed[:3])
            print(
                f"  {name:<22} {'-':>9} {'-':>9} {'-':>13}"
                f"   no result: {len(failed)}/{len(xi)} points fail (xi = {near} A)"
            )
            continue

        transformed = np.asarray(result) / (2 * np.pi) / g_zero_analytic - 1.0
        error = np.max(np.abs(transformed - reference))
        print(
            f"  {name:<22} {call_time * 1e3:6.2f} ms {call_best * 1e3:6.2f} ms"
            f" {call_time / eval_time:12.2f}x {error:11.1e}"
        )


def main():
    model = load_model("sphere")
    for filename, radius, description in DATASETS:
        benchmark_dataset(filename, radius, description, model)


if __name__ == "__main__":
    main()
