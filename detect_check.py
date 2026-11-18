"""Quick check: does Delta-AIC still detect aerosol transmission when the
space-clearance rate zeta is physically realistic?  Single-compartment A-SIR,
synthetic prevalence data (sigma = 0.2%), SIR vs A-SIR fits (eps, zeta fixed)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

BETA, GAMMA, EPS, I0, SIG = 0.4, 0.2, 0.8, 0.01, 0.002
rng = np.random.default_rng(1)

def asir(t, b, g, d, eps, z):
    def f(_, y):
        s, i, w = y
        return [-b*s*i - d*s*w, b*s*i + d*s*w - g*i, eps*(1-w)*i - z*w]
    # seed W on the growth eigenvector (W/I = eps/(r+z), r = b+d*eps/(r+z)-g approx)
    w0 = eps*I0/(BETA-GAMMA+z)
    return solve_ivp(f, (0, t[-1]), [1-I0, I0, w0], t_eval=t, method="LSODA",
                     rtol=1e-8, atol=1e-11).y[1]

def aic(rss, n, k): return n*np.log(rss/n) + 2*k

for z in (0.05, 1.4, 14.0):
    for fa in (0.5, 0.9):
        d_true = BETA*fa*(BETA-GAMMA+z)/EPS          # same early growth as SIR(0.4, 0.2)
        for win in (25, 60):
            t = np.arange(0, win+1, 1.0)
            data = asir(t, BETA*(1-fa), GAMMA, d_true, EPS, z) + rng.normal(0, SIG, t.size)
            best_sir = min((least_squares(lambda p: asir(t, p[0], p[1], 0, EPS, z)-data,
                            x0, bounds=([0,0.01],[3,3])) for x0 in
                            ([0.4,0.2],[0.3,0.15],[0.6,0.3],[0.2,0.1])), key=lambda r: r.cost)
            best_as = min((least_squares(lambda p: asir(t, p[0], p[1], p[2], EPS, z)-data,
                           x0, bounds=([0,0.01,0],[3,3,100])) for x0 in
                           ([0.2,0.2,d_true],[0.4,0.2,0.0],[0.1,0.15,2*d_true])),
                          key=lambda r: r.cost)
            n = t.size
            dA = aic(2*best_sir.cost, n, 2) - aic(2*best_as.cost, n, 3)
            print(f"zeta={z:6.2f}  fa={fa:.1f}  window={win:3d} d  dAIC={dA:8.1f}")
