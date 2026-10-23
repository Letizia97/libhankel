// FFTLog-based Hankel transform.
// Public entry point: sasfit_fftlog(), declared in fft_log.h.
//
// Algorithm: Hamilton (2000), MNRAS 312, 257.  Samples f on a log-spaced
// grid, convolves with the Hankel kernel in DFT space, then returns
// individual query values via a log-log PCHIP spline.  The whole computed
// curve is cached across repeated calls with the same (f, fparams, nu)
// so a q-sweep over many points costs one O(N log N) FFT pass plus cheap
// spline lookups.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <complex>
#include <limits>
#include <vector>

#include "interp_cubic.h"
#include "pocketfft_hdronly.h"
#include "sasfit_lru_cache.h"

// ---------------------------------------------------------------------------
// Complex log-Gamma via the Lanczos approximation (g=7, N=9).
// Coefficients from Numerical Recipes 3e, Table 6.1.
// Accurate to ~15 significant figures for Re(z) >= 0.5, which holds for
// the arguments used in the FFTLog kernel: x = (nu+1)/2 is 0.5 or 1.0.
// ---------------------------------------------------------------------------

namespace {

void clngamma(double re, double im, double *out_re, double *out_im) {
    static const double g[] = {
        0.99999999999980993,  676.5203681218851,  -1259.1392167224028,
        771.32342877765313,  -176.61502916214059,    12.507343278686905,
        -0.13857109526572012,   9.9843695780195716e-6, 1.5056327351493116e-7
    };
    std::complex<double> z(re - 1.0, im);
    std::complex<double> s(g[0], 0.0);
    for (int i = 1; i <= 8; ++i)
        s += g[i] / (z + double(i));
    std::complex<double> t = z + 7.5;
    std::complex<double> r = 0.5 * std::log(2.0 * M_PI)
                           + (z + 0.5) * std::log(t) - t + std::log(s);
    *out_re = r.real();
    *out_im = r.imag();
}

// ---------------------------------------------------------------------------
// Low-ringing kr adjustment (Hamilton 2000, §4).
// Shifts the central k·r product so the m=N/2 kernel coefficient's phase
// lands on a multiple of π/N, suppressing ringing from the DFT's implicit
// circular periodicity.
// ---------------------------------------------------------------------------

double adjusted_kr(int n, double nu, double log_span, double kr0) {
    const double xg = (nu + 1.0) * 0.5;
    const double yg = M_PI * n / (2.0 * log_span);
    double ln_r, phi;
    clngamma(xg, yg, &ln_r, &phi);
    const double a  = std::log(2.0 / kr0) * n / log_span + 2.0 * phi / M_PI;
    const double ia = std::round(a);
    if (a != ia)
        return kr0 * std::exp((a - ia) * log_span / n);
    return kr0;
}

// ---------------------------------------------------------------------------
// Core FFTLog step.
//
// Input:  r_grid (n log-spaced abscissae), a_vals = r * f(r) at those points.
// Output: k_grid (n log-spaced output abscissae), Ak_vals = ∫ a(r) J_nu(k r) k dr.
//
// To recover the standard Hankel transform I(k) = ∫ f(r) J_nu(kr) r dr,
// the caller divides Ak by k.
//
// The output index mapping  out[(n-j) % n]  is required for correct phase
// alignment between the DFT output array and the physical k-grid; using
// the naive out[j] produces an apparent O(1/N) convergence instead of the
// near-spectral convergence FFTLog achieves for smooth inputs.
// ---------------------------------------------------------------------------

struct HankelGrid {
    std::vector<double> k;
    std::vector<double> Ak;
};

HankelGrid fftlog_step(const std::vector<double>& r_grid,
                        const std::vector<double>& a_vals,
                        double nu) {
    const int    n        = static_cast<int>(r_grid.size());
    const int    nc       = n / 2 + 1;
    const double log_span = std::log(r_grid[n-1] / r_grid[0]) * n / (n - 1.0);
    const double kr       = adjusted_kr(n, nu, log_span, 1.0);
    const double k0r0     = kr * std::exp(-log_span);
    const double dfreq    = M_PI / log_span;
    const double phase0   = -2.0 * dfreq * std::log(k0r0 / 2.0);

    // Kernel coefficients u_m (Hamilton 2000, eq. A6)
    const double xg = (nu + 1.0) * 0.5;
    std::vector<std::complex<double>> u(nc);
    for (int m = 0; m < nc; ++m) {
        double ln_r, phi;
        clngamma(xg, m * dfreq, &ln_r, &phi);
        u[m] = std::polar(1.0, m * phase0 + 2.0 * phi);
    }

    // Forward DFT of a(r)
    std::vector<double> buf(a_vals);
    std::vector<std::complex<double>> c(nc);
    pocketfft::r2c({size_t(n)}, {ptrdiff_t(sizeof(double))},
                   {ptrdiff_t(sizeof(std::complex<double>))},
                   0, true, buf.data(), c.data(), 1.0, 1);

    // Multiply by kernel (1/n normalisation absorbed here)
    std::vector<std::complex<double>> cu(nc);
    for (int m = 0; m < nc; ++m)
        cu[m] = c[m] * u[m] / double(n);

    // Inverse DFT
    std::vector<double> out(n);
    pocketfft::c2r({size_t(n)}, {ptrdiff_t(sizeof(std::complex<double>))},
                   {ptrdiff_t(sizeof(double))},
                   0, false, cu.data(), out.data(), 1.0, 1);

    const double k0 = k0r0 / r_grid[0];
    HankelGrid hg;
    hg.k.resize(n);
    hg.Ak.resize(n);
    for (int j = 0; j < n; ++j) {
        hg.k[j]  = k0 * std::exp(j * log_span / n);
        hg.Ak[j] = out[((n - j) % n + n) % n];
    }
    return hg;
}

// ---------------------------------------------------------------------------
// LRU cache of computed curves.
// Each entry stores the full signed Hankel transform values and a log-log
// PCHIP spline of their magnitude, so individual query points are cheap
// spline lookups.  Cache validity is verified via an anchor-ratio check:
// two fixed r-values are re-evaluated on every cache hit and their ratio
// compared against what was recorded at build time.  Using a ratio (not a
// raw value) makes the check invariant to pure scale changes in f, which
// dominate finite-difference gradient evaluations and would otherwise
// cause unnecessary cache misses.
// ---------------------------------------------------------------------------

struct CacheKey {
    double (*func)(double, void*) = nullptr;
    void*  params                 = nullptr;
    double order                  = 0.0;
};

struct CacheEntry {
    cubic_interp_t*     spline    = nullptr;
    std::vector<double> k_pts;
    std::vector<double> ht_vals;  // signed I(k) = Ak/k
    double              anchor_lo = 0.0;
    double              anchor_hi = 0.0;
    double              ratio     = 0.0;

    CacheEntry() = default;

    ~CacheEntry() { drop(); }

    CacheEntry(CacheEntry&& o) noexcept { *this = std::move(o); }
    CacheEntry& operator=(CacheEntry&& o) noexcept {
        if (this != &o) {
            drop();
            spline    = o.spline;  o.spline = nullptr;
            k_pts     = std::move(o.k_pts);
            ht_vals   = std::move(o.ht_vals);
            anchor_lo = o.anchor_lo;
            anchor_hi = o.anchor_hi;
            ratio     = o.ratio;
        }
        return *this;
    }

    CacheEntry(const CacheEntry&)            = delete;
    CacheEntry& operator=(const CacheEntry&) = delete;

    void drop() {
        if (spline) { cubic_interp_destroy(spline); spline = nullptr; }
    }
};

} // namespace

// ---------------------------------------------------------------------------
// Public entry point.
// ---------------------------------------------------------------------------

extern "C" double sasfit_fftlog(double nu, double (*f)(double, void*),
                                 double x, void* fparams) {
    if (!f || !(x > 0.0) || (nu != 0.0 && nu != 1.0))
        return std::numeric_limits<double>::quiet_NaN();

    try {
        static SasfitLruCache<CacheKey, CacheEntry> cache(8);

        CacheKey key{};
        key.func   = f;
        key.params = fparams;
        key.order  = nu;

        // Validate any existing cache entry
        CacheEntry* entry     = cache.lookup(key);
        bool        reuse     = false;
        if (entry && entry->spline && !entry->k_pts.empty()
                  && x >= entry->k_pts.front() && x <= entry->k_pts.back()) {
            const double vlo   = f(entry->anchor_lo, fparams);
            const double vhi   = f(entry->anchor_hi, fparams);
            const double r     = vhi / (vlo + 1e-300);
            const double tol   = 1e-3;
            if (std::fabs(r - entry->ratio)
                    <= tol * (std::fabs(r) + std::fabs(entry->ratio) + 1e-300))
                reuse = true;
        }

        if (!reuse) {
            // Locate the peak of |r·f(r)| and detect compact support.
            // The 121-point coarse probe spans 1e-4 to 1e4 on a log scale.
            std::vector<std::pair<double,double>> probe;
            probe.reserve(121);
            double peak_r = 1.0, peak_v = -1.0;
            for (int i = -60; i <= 60; ++i) {
                const double ri  = std::pow(10.0, i / 15.0);
                const double vi  = std::fabs(ri * f(ri, fparams));
                if (std::isfinite(vi)) {
                    probe.push_back({ri, vi});
                    if (vi > peak_v) { peak_v = vi; peak_r = ri; }
                }
            }

            // Binary-search for the compact-support edge (if any).
            // A value is treated as zero only if it is below 1e-300 and
            // stays there for all larger r — avoids misfires on
            // oscillatory inputs that merely dip low.
            constexpr double ZERO = 1e-300;
            double edge           = -1.0;
            const double srch0    = peak_r * 2.0;
            for (std::size_t i = 0; i < probe.size() && edge < 0.0; ++i) {
                if (probe[i].first < srch0 || probe[i].second >= ZERO) continue;
                bool persist = true;
                for (std::size_t j = i; j < probe.size(); ++j)
                    if (probe[j].second >= ZERO) { persist = false; break; }
                if (persist) {
                    double lo = (i > 0) ? probe[i-1].first : 0.0;
                    double hi = probe[i].first;
                    for (int it = 0; it < 60; ++it) {
                        double mid = 0.5 * (lo + hi);
                        (std::fabs(mid * f(mid, fparams)) < ZERO ? hi : lo) = mid;
                    }
                    edge = hi;
                }
            }

            // Build the log-spaced integration grid.
            // FFTLog needs ~500x the function's natural scale for padding
            // (the periodicity assumption causes large ringing otherwise).
            // The query x also constrains the needed output range via the
            // reciprocal relationship k ~ 1/r.
            const double r_nat_min = peak_r * 1e-6;
            const double r_nat_max = (edge > 0.0 ? edge : peak_r) * 500.0;
            const double r_min     = std::min(r_nat_min, 0.1 / x);
            const double r_max     = std::max(r_nat_max, 10.0  / x);
            const int    n_pts     = 2048;

            std::vector<double> r_grid(n_pts), a_grid(n_pts);
            const double dlnr = std::log(r_max / r_min) / (n_pts - 1);
            for (int j = 0; j < n_pts; ++j) {
                r_grid[j] = r_min * std::exp(j * dlnr);
                const double fv = f(r_grid[j], fparams);
                a_grid[j] = std::isfinite(fv) ? r_grid[j] * fv : 0.0;
            }

            const HankelGrid hg = fftlog_step(r_grid, a_grid, nu);

            // Store signed values and a log-log spline of the magnitude.
            std::vector<double> log_k(n_pts), log_mag(n_pts);
            CacheEntry& ce = cache.insert(key);
            ce.k_pts.resize(n_pts);
            ce.ht_vals.resize(n_pts);
            for (int j = 0; j < n_pts; ++j) {
                const double Ik = hg.Ak[j] / hg.k[j];
                ce.k_pts[j]   = hg.k[j];
                ce.ht_vals[j] = Ik;
                log_k[j]      = std::log(hg.k[j]);
                log_mag[j]    = std::log(std::fabs(Ik) + 1e-300);
            }
            ce.spline = cubic_interp_create(log_k.data(), log_mag.data(), n_pts);

            ce.anchor_lo = 0.5 * peak_r;
            ce.anchor_hi = peak_r;
            ce.ratio     = f(ce.anchor_hi, fparams) / (f(ce.anchor_lo, fparams) + 1e-300);

            entry = &ce;
        }

        if (x < entry->k_pts.front() || x > entry->k_pts.back())
            return 0.0;

        const double log_mag = cubic_interp_eval(entry->spline, std::log(x));

        // Recover sign via binary search on the pre-computed signed values.
        std::size_t lo = 0, hi = entry->k_pts.size() - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (entry->k_pts[mid] < x) lo = mid; else hi = mid;
        }
        const double sign = (entry->ht_vals[lo] < 0.0) ? -1.0 : 1.0;
        return sign * std::exp(log_mag);

    } catch (const std::exception& e) {
        fprintf(stderr, "sasfit_fftlog: %s\n", e.what());
        return std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        fprintf(stderr, "sasfit_fftlog: unknown exception\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
}
