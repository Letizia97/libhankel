"""
The same comparison again, but transforming a table of I(q) instead of a model.

``benchmark.py`` and ``benchmark_real_data.py`` both transform an analytic form
factor, which either route can evaluate at any q it likes. This one transforms
tabulated I(q): a finite list of (q, I) pairs, the shape scattering data
actually arrives in. That changes the comparison completely, because the two
routes need the table in different ways.

* ``SesansTransform`` never asks for a q of its own choosing. ``apply(Iq)`` is
  ``H.T @ Iq``, so it needs I(q) once, on the fixed grid ``q_calc``. Getting a
  table onto that grid is a single vectorised interpolation, and the transform
  itself is unchanged.
* LibHankel chooses its own abscissae, so a table has to be wrapped in a
  callable it can evaluate at arbitrary q. Through the Python interface that
  costs one ``PyObject_CallObject`` per evaluation (``src/py_interface.c``),
  and the strategy decides how many of those there are.

Neither route can invent I(q) outside the table, so both inherit the same
truncation error. That error, not the transform, is what dominates here, and
the third section below measures it on its own.

A caveat on the q range. A SESANS measurement on micron-scale structure needs
I(q) down to ~1e-5 1/A, which is below what a SANS instrument reaches: the D22
curves distributed with sasmodels start at 0.04 1/A. So the main table here is
synthetic and deliberately generous (1e-06 to 1e+01 1/A), to make the transform
error visible rather than drowning it in truncation. Section three then shows
what a real instrument window does to the answer.

The spin-echo lengths are the real ones from ``sphere_isis.ses``.

Requires ``sasmodels`` and ``sasdata``:

    python -m pip install sasmodels sasdata
    python examples/sesans/benchmark_tabulated.py
"""

import bisect
import math
import time

import libhankel
import numpy as np
from sasdata import data_path
from sasmodels.core import load_model
from sasmodels.data import load_data
from sasmodels.direct_model import DirectModel

SES_FILE = data_path / "sesans_data" / "sphere_isis.ses"

RADIUS = 1000.0  # A
CONTRAST = 1.0

# The reference table: log-spaced, wider than any real instrument.
Q_MIN, Q_MAX = 1e-6, 1e1
N_TABLE = 1000

# Fewer repeats than the other benchmarks: the Python callback makes some of
# these strategies hundreds of times slower, and the timings are stable anyway.
REPEATS = 5

# G(0) for a sphere, in the 1 / 2 pi convention used by SesansTransform.
G_ZERO = 2 * np.pi * CONTRAST**2 * RADIUS**4

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

# Table sizes for the density sweep, and q windows for the truncation study.
# The instrument windows are the q ranges of real configurations: the first is
# the D22 curve shipped with sasmodels, the last is a generous SANS + USANS
# merge. None of them reaches the 1e-5 that this sample needs.
TABLE_SIZES = [100, 300, 1000, 3000, 10000]
INSTRUMENT_WINDOWS = [
    ("D22 as shipped", 0.0399, 0.6128),
    ("typical SANS", 3e-3, 0.5),
    ("SANS + USANS", 1e-4, 0.5),
    ("USANS to low q", 1e-5, 0.5),
    ("generous synthetic", Q_MIN, Q_MAX),
]


def sphere_iq(q, radius=RADIUS, contrast=CONTRAST):
    """
    The same I(q) as LibHankel's built-in "sphere", so that any difference
    below comes from the tabulation and not from a different model.
    """
    return (
        contrast * 4.0 * np.pi * (np.sin(q * radius) - q * radius * np.cos(q * radius)) / q**3
    ) ** 2


def analytic_sesans_sphere(xi, radius=RADIUS):
    """P(xi) / G(0) for a homogeneous sphere, the reference for both routes."""
    u = np.clip(xi / (2.0 * radius), 0.0, 1.0)
    root = np.sqrt(1.0 - u**2)
    log_term = np.log(np.maximum(u / (1.0 + root), 1e-15))
    return (root * (2.0 + u**2) + u**2 * (4.0 - u**2) * log_term - 2.0) / 2.0


class Table:
    """
    Tabulated I(q) with log-log interpolation, the usual choice for SAS data.

    Outside the table I(q) is zero: a table cannot say anything about q it does
    not cover, and pinning the end values instead would put a non-integrable
    tail into the transform. The call counters are what make the cost of the
    LibHankel route legible.
    """

    def __init__(self, q_min=Q_MIN, q_max=Q_MAX, n=N_TABLE):
        self.q = np.logspace(np.log10(q_min), np.log10(q_max), n)
        self.iq = sphere_iq(self.q)
        # Pure Python lists and bisect, rather than np.interp on a scalar:
        # numpy's scalar path would add overhead that has nothing to do with
        # the transform, and a real binding would not pay it.
        self.log_q = [math.log(value) for value in self.q]
        self.log_iq = [math.log(value) if value > 0.0 else -700.0 for value in self.iq]
        self.q_lo, self.q_hi = self.q[0], self.q[-1]
        self.calls = 0
        self.outside = 0

    def __call__(self, q, params):
        """Scalar interpolation, called once per abscissa by LibHankel."""
        self.calls += 1
        if q <= self.q_lo or q >= self.q_hi:
            self.outside += 1
            return 0.0
        i = bisect.bisect_right(self.log_q, math.log(q)) - 1
        lo, hi = self.log_q[i], self.log_q[i + 1]
        weight = (math.log(q) - lo) / (hi - lo)
        return math.exp(self.log_iq[i] * (1.0 - weight) + self.log_iq[i + 1] * weight)

    def on_grid(self, q_grid):
        """Vectorised interpolation onto a fixed grid, for the sasmodels route."""
        values = np.exp(np.interp(np.log(q_grid), self.log_q, self.log_iq))
        return np.where((q_grid > self.q_lo) & (q_grid < self.q_hi), values, 0.0)

    def g_zero(self):
        """G(0) from the table alone: integral of q I(q) / 2 pi, in log q."""
        return np.trapezoid(self.q**2 * self.iq, np.log(self.q)) / (2 * np.pi)

    def reset(self):
        self.calls = self.outside = 0


def timings(call, repeats=REPEATS):
    """Median and best wall-clock time over *repeats* calls, plus the result."""
    measured = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = call()
        measured.append(time.perf_counter() - start)
    return float(np.median(measured)), min(measured), result


def sasmodels_route(transform, table):
    """Interpolate the table onto q_calc, then apply H. One model-free evaluation."""
    iq = table.on_grid(transform.q_calc)
    return transform.apply(iq), np.dot(transform._H0, iq)


def main_comparison(xi, transform, reference):
    table = Table()
    print(f"\ntabulated I(q): {N_TABLE} points, {Q_MIN:g} to {Q_MAX:g} 1/A, log-log interpolated")
    print(
        f"  G(0) from the table is off by {abs(table.g_zero() / G_ZERO - 1):.1e}"
        f" relative to the analytic 2 pi eta^2 R^4"
    )

    eval_time, eval_best, (polarisation, g_zero) = timings(
        lambda: sasmodels_route(transform, table)
    )
    sasmodels_error = np.max(np.abs(polarisation / g_zero - reference))
    print(
        f"  sasmodels  interpolate + apply    {eval_time * 1e3:7.2f} ms"
        f"  (best {eval_best * 1e3:.2f})   max error {sasmodels_error:.1e}"
    )
    print(
        f"  {'libhankel strategy':<22} {'median':>9} {'best':>9} {'vs sasmodels':>13}"
        f" {'max error':>11} {'I(q) calls':>12} {'outside':>9}"
    )

    for name, strategy_params in STRATEGIES:
        table.reset()
        try:
            call_time, call_best, result = timings(
                lambda: libhankel.hankel_transform(
                    0,
                    table,
                    xi,
                    np.array([RADIUS, CONTRAST]),
                    name,  # noqa: B023
                    strategy_params,  # noqa: B023
                )
            )
        except RuntimeError:
            # The counters still hold whatever the strategy managed before it
            # gave up, which is the interesting part of a failure here.
            print(
                f"  {name:<22} {'-':>9} {'-':>9} {'-':>13} {'no result':>11}"
                f" {table.calls:12d} {100.0 * table.outside / table.calls:8.0f}%"
            )
            continue

        transformed = np.asarray(result) / (2 * np.pi) / table.g_zero() - 1.0
        error = np.max(np.abs(transformed - reference))
        per_call = table.calls // REPEATS
        outside = 100.0 * table.outside / table.calls
        print(
            f"  {name:<22} {call_time * 1e3:6.2f} ms {call_best * 1e3:6.2f} ms"
            f" {call_time / eval_time:12.2f}x {error:11.1e} {per_call:12d} {outside:8.0f}%"
        )


def density_sweep(xi, model, reference):
    print("\nhow many tabulated points are needed, over the full 1e-06 to 1e+01 range")
    print(f"  {'points':>8} {'sasmodels':>12} {'DHT_Key_201':>13} {'DHT_Anderson_801':>18}")
    for n in TABLE_SIZES:
        table = Table(n=n)
        transform = DirectModel(load_data(str(SES_FILE)), model).resolution
        polarisation, g_zero = sasmodels_route(transform, table)
        errors = [np.max(np.abs(polarisation / g_zero - reference))]
        for name in ("DHT_Key_201", "DHT_Anderson_801"):
            result = libhankel.hankel_transform(
                0, table, xi, np.array([RADIUS, CONTRAST]), name, {}
            )
            transformed = np.asarray(result) / (2 * np.pi) / table.g_zero() - 1.0
            errors.append(np.max(np.abs(transformed - reference)))
        print(f"  {n:8d} {errors[0]:12.1e} {errors[1]:13.1e} {errors[2]:18.1e}")


def truncation_study(xi, model, reference):
    print("\nwhat the measured q window costs, at a fixed 1000 points per window")
    print(f"  {'window':<22} {'q range (1/A)':>22} {'sasmodels':>12} {'DHT_Key_201':>13}")
    for label, q_min, q_max in INSTRUMENT_WINDOWS:
        table = Table(q_min=q_min, q_max=q_max, n=1000)
        transform = DirectModel(load_data(str(SES_FILE)), model).resolution
        polarisation, g_zero = sasmodels_route(transform, table)
        sasmodels_error = np.max(np.abs(polarisation / g_zero - reference))
        result = libhankel.hankel_transform(
            0, table, xi, np.array([RADIUS, CONTRAST]), "DHT_Key_201", {}
        )
        transformed = np.asarray(result) / (2 * np.pi) / table.g_zero() - 1.0
        error = np.max(np.abs(transformed - reference))
        print(f"  {label:<22} {q_min:9.1e} to {q_max:8.1e} {sasmodels_error:12.1e} {error:13.1e}")


def main():
    data = load_data(str(SES_FILE))
    xi = data.x
    reference = analytic_sesans_sphere(xi)
    model = load_model("sphere")
    transform = DirectModel(data, model).resolution

    print(f"sphere of radius {RADIUS:g} A, spin-echo lengths from {SES_FILE.name}")
    print(f"  {len(xi)} values, {xi.min():.0f} to {xi.max():.0f} A")
    print(
        f"  SesansTransform q_calc: {len(transform.q_calc)} points,"
        f" {transform.q_calc[0]:.2e} to {transform.q_calc[-1]:.2e} 1/A"
    )

    main_comparison(xi, transform, reference)
    density_sweep(xi, model, reference)
    truncation_study(xi, model, reference)


if __name__ == "__main__":
    main()
