"""
Swap the sasmodels SESANS transform for a LibHankel filter and measure what
changes.

``benchmark_real_data.py`` compares ``libhankel.hankel_transform`` against the
sasmodels route as two separate programs computing the same integral. It has to
supply G(0) analytically, and it never applies the analyser acceptance mask, so
the two sides are not quite computing the same thing.

This benchmark closes both gaps by making the filter a sasmodels resolution
object instead of a separate program. ``DHTSesansTransform`` in ``dht_sesans``
exposes the ``q_calc`` and ``apply`` that ``DirectModel._calc_theory`` uses, so
the fit driver, the parameter handling, the unit conversions, the acceptance
masking and the polydispersity all stay exactly as they are, and only the
quadrature underneath changes::

    original = use_dht_transform("DHT_Key_201")
    calculator = DirectModel(data, model)   # now uses the filter
    calculator(radius=1000.0, ...)          # returns P(xi) as before

Three things are reported per dataset:

* build time, which is where the dense-matrix route spends most of its setup:
  it allocates and fills an (n_q, n_xi) matrix of Bessel values;
* evaluation time, the cost per fit iteration, which is what the thousands of
  iterations in a fit actually multiply;
* accuracy against the analytic sphere, so a speedup that came from doing less
  arithmetic is distinguishable from one that came from doing it worse.

Requires ``sasmodels`` and ``sasdata``:

    python -m pip install sasmodels sasdata
    python examples/sesans/benchmark_dht_transform.py
"""

import time

import numpy as np
from dht_sesans import FILTERS, DHTSesansTransform
from sasdata import data_path
from sasmodels.core import load_model
from sasmodels.data import load_data
from sasmodels.direct_model import DirectModel, call_kernel
from sasmodels.sesans import SesansTransform

SESANS_DIR = data_path / "sesans_data"

DATASETS = [
    ("sphere_isis.ses", 1000.0, "PMMA in deuterated decalin, ISIS, time-of-flight"),
    ("sphere2micron.ses", 10000.0, "2 um polystyrene, Delft, monochromatic"),
    ("sphere2micron_long.ses", 10000.0, "the same sample, out to 20 um"),
]

MODEL_PARS = dict(sld=1.41, sld_solvent=2.70, scale=1.0, background=0.0)

REPEATS = 15


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


def sesans_arguments(data):
    """
    The five arguments ``_make_sesans_transform`` derives from a data file.

    Both transforms take these, in this order. Reproduced here rather than
    called through sasmodels so that the transform can be built on its own,
    away from a DirectModel, for the build-time measurement.
    """
    theta_max, _ = data.sample.zacceptance
    wavelength = np.asarray(data.source.wavelength, dtype=float)
    zaccept = 2 * np.pi / np.max(wavelength) * np.sin(theta_max)
    return data.x, np.asarray(data.x, dtype=float), wavelength, zaccept, 10000000


def evaluate(data, model, transform, radius):
    """
    Run a model through *transform* and return (median time, best, P/G(0)).

    ``DirectModel`` builds its own transform in the constructor; assigning over
    it and clearing the cached kernel is what makes the swap take effect,
    because ``_calc_theory`` reads ``resolution.q_calc`` only when it has no
    kernel yet.
    """
    calculator = DirectModel(data, model)
    if transform is not None:
        calculator.resolution = transform
        calculator._kernel = None
    calculator(radius=radius, **MODEL_PARS)  # warm up: the first call compiles the kernel

    median, best, polarisation = timings(lambda: calculator(radius=radius, **MODEL_PARS))
    g_zero = np.sum(calculator.resolution._H0 * calculator.Iq_calc)
    return median, best, polarisation / g_zero, calculator.resolution


def benchmark_dataset(filename, radius, description, model):
    data = load_data(str(SESANS_DIR / filename))
    xi = np.asarray(data.x, dtype=float)
    reference = analytic_sesans_sphere(xi, radius)
    args = sesans_arguments(data)

    print(f"\n{filename}   {description}")
    print(
        f"  {len(xi)} spin-echo lengths, {xi.min():.0f} to {xi.max():.0f} A;  radius {radius:g} A"
    )
    print(
        f"  {'transform':<20} {'build':>9} {'per eval':>9} {'best':>9}"
        f" {'q points':>9} {'speedup':>8} {'max error':>11}"
    )

    build, _, _ = timings(lambda: SesansTransform(*args), repeats=3)
    median, best, result, resolution = evaluate(data, model, None, radius)
    error = np.max(np.abs(result - reference))
    print(
        f"  {'sasmodels dense':<20} {build * 1e3:6.2f} ms {median * 1e3:6.2f} ms"
        f" {best * 1e3:6.2f} ms {len(resolution.q_calc):9d} {1.0:7.1f}x {error:11.1e}"
    )
    baseline = median

    for strategy in FILTERS:
        build, _, _ = timings(
            lambda: DHTSesansTransform(*args, strategy=strategy),  # noqa: B023
            repeats=3,
        )
        transform = DHTSesansTransform(*args, strategy=strategy)
        median, best, result, resolution = evaluate(data, model, transform, radius)
        error = np.max(np.abs(result - reference))
        print(
            f"  {strategy:<20} {build * 1e3:6.2f} ms {median * 1e3:6.2f} ms"
            f" {best * 1e3:6.2f} ms {len(resolution.q_calc):9d}"
            f" {baseline / median:7.1f}x {error:11.1e}"
        )


def check_g0(model):
    """
    How well each route computes G(0), the one thing a filter cannot do directly.

    ``apply`` returns ``G(xi) - G(0)``, and G(0) is the transform at xi = 0. A
    filter evaluates at ``node/xi``, so it cannot go there at all;
    ``SesansTransform`` gets G(0) free from its dense grid. ``DHTSesansTransform``
    instead runs a trapezoid rule in log q over the union of the abscissae it is
    already using, which costs no extra model evaluations. This checks that the
    substitute is good enough, against a reference from a deliberately
    over-resolved grid.
    """
    data = load_data(str(SESANS_DIR / "sphere_isis.ses"))
    args = sesans_arguments(data)
    pars = dict(radius=1000.0, **MODEL_PARS)
    pars["background"] = 0.0

    def g0_on(q):
        Iq = call_kernel(model.make_kernel((q,)), pars)
        return Iq, np.trapezoid(q**2 * Iq, np.log(q)) / (2 * np.pi)

    fine = np.exp(np.linspace(np.log(1e-9), np.log(1e4), 400000))
    _, reference = g0_on(fine)

    print("\nG(0), the integral a filter cannot reach directly")
    print(f"  {'route':<20} {'q points':>9} {'relative error':>15}")
    print(f"  {'reference':<20} {len(fine):9d}   400k point log grid")

    transform = SesansTransform(*args)
    Iq = call_kernel(model.make_kernel((transform.q_calc,)), pars)
    error = np.sum(transform._H0 * Iq) / reference - 1.0
    print(f"  {'sasmodels dense':<20} {len(transform.q_calc):9d} {error:15.1e}")

    for strategy in FILTERS:
        transform = DHTSesansTransform(*args, strategy=strategy)
        Iq = call_kernel(model.make_kernel((transform.q_calc,)), pars)
        error = np.sum(transform._H0 * Iq) / reference - 1.0
        print(f"  {strategy:<20} {len(transform.q_calc):9d} {error:15.1e}")


def check_drop_in(model):
    """
    Show that the swap needs no changes to the calling code.

    ``use_dht_transform`` replaces the one function that decides which
    transform a SESANS DirectModel gets. Everything downstream -- the
    constructor, the call, the returned P(xi) -- is untouched.
    """
    from dht_sesans import use_dht_transform
    from sasmodels import direct_model

    data = load_data(str(SESANS_DIR / "sphere_isis.ses"))
    plain = DirectModel(data, model)
    dense = plain(radius=1000.0, **MODEL_PARS)

    original = use_dht_transform("DHT_Key_201")
    try:
        calculator = DirectModel(data, model)
        filtered = calculator(radius=1000.0, **MODEL_PARS)
        filtered_g0 = np.sum(calculator.resolution._H0 * calculator.Iq_calc)
    finally:
        direct_model._make_sesans_transform = original

    reference = analytic_sesans_sphere(np.asarray(data.x, dtype=float), 1000.0)
    dense_g0 = np.sum(plain.resolution._H0 * plain.Iq_calc)
    dense_error = np.max(np.abs(dense / dense_g0 - reference))
    filtered_error = np.max(np.abs(filtered / filtered_g0 - reference))

    print("\ndrop-in check: DirectModel patched to use DHT_Key_201")
    print(f"  transform in use          {type(calculator.resolution).__name__}")
    print(f"  P(xi) returned            shape {filtered.shape}, as before")
    print(f"  max |P/G0 - analytic|     {dense_error:.1e} dense")
    print(f"                            {filtered_error:.1e} filter")
    print("  calling code unchanged:   DirectModel(data, model)(radius=..., **pars)")


def main():
    model = load_model("sphere")
    print("Replacing the sasmodels SESANS transform with a LibHankel digital filter.")
    print("Speedup is on the per-evaluation time, which is what a fit multiplies.")
    for filename, radius, description in DATASETS:
        benchmark_dataset(filename, radius, description, model)
    check_g0(model)
    check_drop_in(model)


if __name__ == "__main__":
    main()
