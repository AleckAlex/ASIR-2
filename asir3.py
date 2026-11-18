"""
asir3.py  --  SIR-1, SIR-3 and A-SIR-3 models with calibration to a common
              early growth rate, plus the mitigation grid of the paper.

Author : A.H. Alexopoulos (CERTH), with Python drafting help
Python : 3.10+  (numpy, scipy, matplotlib)

Notes for a Fortran programmer
------------------------------
* Arrays are 0-based:  y[0] is the first element (Fortran y(1)).
* numpy arrays behave like Fortran arrays for arithmetic: a*b is element-wise.
* A "dict" (e.g. MITIGATIONS below) is a lookup table keyed by a name.
* solve_ivp(f, ...) plays the role of a call to an ODE library (e.g. DASSL /
  LSODE); f(t, y) is the RHS subroutine and must return dy/dt as an array.
* Run:  python asir3.py      -> prints tables, writes PNG figures and a CSV.

Model hierarchy (see manuscript)
--------------------------------
SIR-1   : classical single-compartment SIR with effective beta_eff, gamma.
SIR-3   : three spaces (0 outdoors, 1 home, 2 work), direct transmission only.
          With fast exchange it collapses to SIR-1  (checked in fig A).
A-SIR-3 : SIR-3 plus infected indoor spaces W1, W2 (aerosol pathway).

Calibration
-----------
For aerosol fraction fa the direct rate is beta_eff*(1-fa) and delta is found
by root-finding so that the early exponential growth rate r of A-SIR-3 equals
that of the SIR-1 reference, r = beta_eff - gamma.  At the exponential phase
the aerosol share of the force of infection is then exactly fa.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 1. Global settings
# ---------------------------------------------------------------------------
MIXING   = "frequency"   # "frequency": beta*S_i*I_i/N_i   (recommended)
                         # "density"  : beta*S_i*I_i       (as in draft v3)
KAPPA    = 0.2           # outdoor/indoor direct-transmission ratio (beta0 = KAPPA*beta_in)
GAMMA    = 0.1           # recovery rate, 1/day  (R0 = 3, 5 for beta_eff = 0.3, 0.5)
EPS      = (0.5, 0.8)    # aerosol emission rate, home / work  (1/day)
ZETA     = (0.05, 0.05)  # space clearance rate, home / work  (1/day)
I_INIT   = 0.01          # initial infected fraction
T_END    = 360.0         # days
F_BASE   = np.array([2.0, 12.0, 10.0]) / 24.0   # time budget: outdoors, home, work

AEROSOL_FRACTIONS = [0.0, 0.2, 0.5, 0.9]
BETA_EFF_LIST     = [0.3, 0.5]

# Mitigation multipliers (Table 1).  Work-from-home is now implemented through
# the time budget (23 h home, 1 h outdoors, 0 h work) together with beta x 0.2.
# NO epsilon x 1.92: the longer home occupancy raises I_1 and hence home
# contamination automatically (I_1/I goes from 12/24 to 23/24 = 1.92x).
MITIGATIONS = {
    #  name          beta   delta  eps   zeta   time budget (h: out, home, work)
    "None":        (1.0,  1.0,  1.0,  1.0,  None),
    "Barriers":    (0.4,  1.0,  1.0,  1.0,  None),
    "Distancing":  (0.5,  1.0,  1.0,  1.0,  None),
    "Ventilation": (1.0,  1.0,  1.0, 10.0,  None),
    "Masking":     (0.1,  0.4,  0.6,  1.0,  None),
    "WFH":         (0.2,  1.0,  1.0,  1.0,  (1.0, 23.0, 0.0)),
}

# State vector layout (index -> variable); Fortran-style named indices
iS0, iS1, iS2, iI0, iI1, iI2, iR0, iR1, iR2, iW1, iW2, iCd, iCa = range(13)
NSTATE = 13          # Cd, Ca = cumulative infections via direct / aerosol path


# ---------------------------------------------------------------------------
# 2. Parameter set
# ---------------------------------------------------------------------------
def exchange_rates(f):
    """q_in[i], q_out[i] for i = 1, 2 from the time budget f = (f0, f1, f2).
    One visit per day: residence time tau_i = f_i (days), q_out = 1/tau_i,
    steady state q_in*S0 = q_out*S_i  ->  q_in = (f_i/f0) q_out = 1/f0."""
    q_in, q_out = np.zeros(3), np.zeros(3)
    for i in (1, 2):
        if f[i] > 0.0:
            q_out[i] = 1.0 / f[i]
            q_in[i]  = 1.0 / f[0]
        else:                       # space not used (e.g. work under WFH)
            q_out[i] = 24.0         # empties within ~1 h
            q_in[i]  = 0.0
    return q_in, q_out


def make_params(beta_dir_eff, delta, f=F_BASE, mult=(1, 1, 1, 1)):
    """Build the parameter dictionary.
    beta_dir_eff : effective direct rate seen by SIR-1 (before mitigation)
    delta        : aerosol infection rate (same for home and work)
    mult         : (beta, delta, eps, zeta) multipliers of a mitigation."""
    mb, md, me, mz = mult
    # beta_eff = sum_i f_i beta_i (frequency mixing) with beta_0 = KAPPA*beta_in
    if MIXING == "frequency":
        w = F_BASE[1] + F_BASE[2] + KAPPA * F_BASE[0]
    else:  # density mixing: beta_eff = sum_i f_i^2 beta_i
        w = F_BASE[1]**2 + F_BASE[2]**2 + KAPPA * F_BASE[0]**2
    beta_in = beta_dir_eff / w
    q_in, q_out = exchange_rates(f)
    return dict(beta=np.array([KAPPA * beta_in, beta_in, beta_in]) * mb,
                delta=np.array([0.0, delta, delta]) * md,
                eps=np.array([0.0, EPS[0], EPS[1]]) * me,
                zeta=np.array([1.0, ZETA[0], ZETA[1]]) * mz,
                gamma=GAMMA, q_in=q_in, q_out=q_out, f=np.asarray(f))


# ---------------------------------------------------------------------------
# 3. Right-hand side of A-SIR-3   (SIR-3 is the special case delta = 0)
# ---------------------------------------------------------------------------
def rhs(t, y, p):
    S = y[iS0:iS2 + 1]; I = y[iI0:iI2 + 1]; R = y[iR0:iR2 + 1]
    W = np.array([0.0, y[iW1], y[iW2]])            # no W outdoors
    N = np.maximum(S + I + R, 1e-12)               # occupancy of each space

    if MIXING == "frequency":
        inf_dir = p["beta"] * S * I / N
    else:
        inf_dir = p["beta"] * S * I
    inf_aer = p["delta"] * S * W
    rec     = p["gamma"] * I

    dy = np.zeros(NSTATE)
    for X, iX in ((S, iS0), (I, iI0), (R, iR0)):   # human exchange flows
        for i in (1, 2):
            flow = p["q_in"][i] * X[0] - p["q_out"][i] * X[i]   # outdoors -> i
            dy[iX + i] += flow
            dy[iX]     -= flow
    dy[iS0:iS2 + 1] += -inf_dir - inf_aer
    dy[iI0:iI2 + 1] +=  inf_dir + inf_aer - rec
    dy[iR0:iR2 + 1] +=  rec
    for i, iW in ((1, iW1), (2, iW2)):              # infected spaces
        dy[iW] = p["eps"][i] * (1.0 - W[i]) * I[i] - p["zeta"][i] * W[i]
    dy[iCd] = inf_dir.sum()
    dy[iCa] = inf_aer.sum()
    return dy


def sir1(beta, gamma=GAMMA, i0=I_INIT, t_end=T_END):
    """Classical single-compartment SIR reference."""
    def f(t, y):
        s, i = y
        return [-beta * s * i, beta * s * i - gamma * i]
    t = np.linspace(0, t_end, int(t_end * 4) + 1)
    sol = solve_ivp(f, (0, t_end), [1 - i0, i0], t_eval=t, method="Radau",
                    rtol=1e-8, atol=1e-11)
    return sol.t, sol.y[1]


# ---------------------------------------------------------------------------
# 4. Linearisation at the disease-free state: growth rate, R0, eigenvector
# ---------------------------------------------------------------------------
def linear_matrices(p):
    """Return F (new infections) and V (transitions) for x = [I0,I1,I2,W1,W2]
    at the disease-free state (S_i = N_i = f_i, Z = 1), so dx/dt = (F - V) x."""
    f, b, d, e, z = p["f"], p["beta"], p["delta"], p["eps"], p["zeta"]
    qi, qo, g = p["q_in"], p["q_out"], p["gamma"]
    bb = b if MIXING == "frequency" else b * f
    F = np.zeros((5, 5)); V = np.zeros((5, 5))
    for i in range(3):
        F[i, i] = bb[i]
        V[i, i] += g
    F[1, 3] = d[1] * f[1]; F[2, 4] = d[2] * f[2]          # aerosol infections
    V[0, 0] += qi[1] + qi[2]; V[0, 1] -= qo[1]; V[0, 2] -= qo[2]
    V[1, 1] += qo[1]; V[1, 0] -= qi[1]
    V[2, 2] += qo[2]; V[2, 0] -= qi[2]
    V[3, 3] = z[1]; V[3, 1] = -e[1]                        # contamination of spaces
    V[4, 4] = z[2]; V[4, 2] = -e[2]
    return F, V


def growth_rate(p):
    F, V = linear_matrices(p)
    lam, vec = np.linalg.eig(F - V)
    k = np.argmax(lam.real)
    v = np.abs(vec[:, k].real)
    return lam[k].real, v


def R0_of(p):
    F, V = linear_matrices(p)
    return max(abs(np.linalg.eigvals(F @ np.linalg.inv(V))))


def calibrate_delta(beta_eff, fa):
    """delta such that A-SIR-3 grows at r = beta_eff - gamma (SIR-1 rate)."""
    if fa == 0.0:
        return 0.0
    r_target = beta_eff - GAMMA
    func = lambda d: growth_rate(make_params(beta_eff * (1 - fa), d))[0] - r_target
    hi = 1.0
    while func(hi) < 0:
        hi *= 2.0
    return brentq(func, 0.0, hi, xtol=1e-12)


# ---------------------------------------------------------------------------
# 5. Simulation and metrics
# ---------------------------------------------------------------------------
def initial_state(p_seed):
    """Seed along the dominant growth eigenvector of the unmitigated system,
    scaled to total infected I_INIT; humans distributed by the baseline budget."""
    _, v = growth_rate(p_seed)
    v = v / v[:3].sum() * I_INIT
    y0 = np.zeros(NSTATE)
    y0[iI0:iI2 + 1] = v[:3]
    y0[iW1], y0[iW2] = v[3], v[4]
    y0[iS0:iS2 + 1] = F_BASE - v[:3]
    return y0


def simulate(p, y0, t_end=T_END):
    t = np.linspace(0, t_end, int(t_end * 4) + 1)
    sol = solve_ivp(rhs, (0, t_end), y0, args=(p,), t_eval=t, method="Radau",
                    max_step=0.5, rtol=1e-8, atol=1e-11)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def metrics(t, y):
    I = y[iI0:iI2 + 1].sum(axis=0)
    peak = I.max()
    auc90 = np.trapezoid(I[t <= 90], t[t <= 90])
    # duration: time I stays above the initial prevalence (after first exceeding it)
    above = I > I_INIT * 1.0001
    dur = t[above][-1] - t[above][0] if above.any() else 0.0
    cd, ca = y[iCd, -1], y[iCa, -1]
    aer_share = ca / (cd + ca) if cd + ca > 0 else 0.0
    return dict(peak=100 * peak, t_peak=t[np.argmax(I)], auc90=auc90,
                duration=dur, attack=100 * (cd + ca), aer_share=100 * aer_share)


# ---------------------------------------------------------------------------
# 6. Main
# ---------------------------------------------------------------------------
def main():
    import argparse
    global ZETA
    ap = argparse.ArgumentParser(description="SIR-1 / SIR-3 / A-SIR-3 runs")
    ap.add_argument("--zeta-home", type=float, default=ZETA[0], help="1/day")
    ap.add_argument("--zeta-work", type=float, default=ZETA[1], help="1/day")
    ap.add_argument("--tag", default="", help="suffix for output file names")
    args = ap.parse_args()
    ZETA = (args.zeta_home, args.zeta_work)
    tag = ("_" + args.tag) if args.tag else ""
    print(f"MIXING = {MIXING},  KAPPA = {KAPPA},  gamma = {GAMMA},  "
          f"zeta = {ZETA} 1/day\n")

    # --- (a) calibration table ------------------------------------------------
    print("Calibration (delta found so that A-SIR-3 early growth = SIR-1 growth)")
    print(f"{'beta_eff':>8} {'fa':>5} {'delta':>9} {'r':>8} {'r_SIR1':>8} {'R0':>7}")
    cal = {}
    for be in BETA_EFF_LIST:
        for fa in AEROSOL_FRACTIONS:
            d = calibrate_delta(be, fa)
            p = make_params(be * (1 - fa), d)
            r, _ = growth_rate(p)
            cal[(be, fa)] = d
            print(f"{be:8.2f} {fa:5.2f} {d:9.4f} {r:8.4f} {be - GAMMA:8.4f} {R0_of(p):7.2f}")

    # --- (b) SIR-3 vs SIR-1 check --------------------------------------------
    figA, axA = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    print("\nSIR-3 vs SIR-1 (delta = 0):")
    for k, be in enumerate(BETA_EFF_LIST):
        p = make_params(be, 0.0)
        t, y = simulate(p, initial_state(p))
        I3 = y[iI0:iI2 + 1].sum(axis=0)
        t1, I1 = sir1(be)
        err = np.max(np.abs(I3 - I1))
        print(f"  beta_eff = {be}: peak SIR-3 = {100*I3.max():.2f}%, "
              f"SIR-1 = {100*I1.max():.2f}%, max |dI| = {100*err:.3f} pp")
        axA[k].plot(t1, 100 * I1, color="#2a78d6", lw=2, label="SIR-1")
        axA[k].plot(t, 100 * I3, color="#eb6834", lw=2, ls="--", label="SIR-3")
        axA[k].set_title(f"β_eff = {be}, R₀ = {be/GAMMA:.0f}", fontsize=10)
        axA[k].set_xlim(0, 200); axA[k].set_xlabel("Days")
    axA[0].set_ylabel("Infected (%)"); axA[0].legend(frameon=False)
    for a in axA:
        a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25)
    figA.tight_layout(); figA.savefig(f"figA_sir3_vs_sir1{tag}.png", dpi=200)

    # --- (c) mitigation grid ---------------------------------------------------
    styles = {"None": ("#55554f", "-"), "Barriers": ("#2a78d6", "--"),
              "Distancing": ("#1baf7a", ":"), "Ventilation": ("#eb6834", "-."),
              "Masking": ("#e87ba4", (0, (1, 1))), "WFH": ("#eda100", "-")}
    fas = [0.0, 0.5, 0.9]
    figB, axB = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True, sharey=True)
    rows = []
    for r_, be in enumerate(BETA_EFF_LIST):
        for c_, fa in enumerate(AEROSOL_FRACTIONS):
            d = cal[(be, fa)]
            p_seed = make_params(be * (1 - fa), d)
            y0 = initial_state(p_seed)
            for name, (mb, md, me, mz, budget) in MITIGATIONS.items():
                f = F_BASE if budget is None else np.array(budget) / 24.0
                p = make_params(be * (1 - fa), d, f=f, mult=(mb, md, me, mz))
                t, y = simulate(p, y0)
                m = metrics(t, y)
                rows.append((be, fa, name, m))
                if fa in fas:
                    col, ls = styles[name]
                    ax = axB[r_, fas.index(fa)]
                    ax.plot(t, 100 * y[iI0:iI2 + 1].sum(axis=0), color=col,
                            ls=ls, lw=2 if name != "WFH" else 1.6, label=name)
    for r_, be in enumerate(BETA_EFF_LIST):
        for c_, fa in enumerate(fas):
            ax = axB[r_, c_]
            ax.set_title(f"β_eff = {be} (R₀ = {be/GAMMA:.0f}) | {int(fa*100)}% aerosol",
                         fontsize=10)
            ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=0.25)
            ax.set_xlim(0, 360)
        axB[r_, 0].set_ylabel("Infected (%)")
    for ax in axB[1]:
        ax.set_xlabel("Days")
    h, l = axB[0, 0].get_legend_handles_labels()
    figB.legend(h, l, loc="upper center", ncol=6, frameon=False)
    figB.suptitle(f"ζ_home = {ZETA[0]:g}, ζ_work = {ZETA[1]:g} day⁻¹", y=0.905, fontsize=9, color="#55554f")
    figB.tight_layout(rect=(0, 0, 1, 0.94))
    figB.savefig(f"figB_mitigation_grid{tag}.png", dpi=200)

    # --- (d) metrics table ------------------------------------------------------
    print("\nMitigation metrics")
    hdr = (f"{'beta':>5} {'fa':>4} {'measure':<12} {'peak%':>7} {'red.%':>6} "
           f"{'t_pk':>6} {'AUC90':>7} {'dur':>6} {'attack%':>8} {'aer.inf%':>9}")
    print(hdr)
    with open(f"mitigation_metrics{tag}.csv", "w") as fh:
        fh.write("beta_eff,aerosol_fraction,measure,peak_pct,peak_reduction_pct,"
                 "t_peak,auc90,duration_days,attack_pct,aerosol_share_of_infections_pct\n")
        base = {}
        for be, fa, name, m in rows:
            if name == "None":
                base[(be, fa)] = m["peak"]
            red = 100 * (1 - m["peak"] / base[(be, fa)])
            print(f"{be:5.2f} {fa:4.1f} {name:<12} {m['peak']:7.2f} {red:6.1f} "
                  f"{m['t_peak']:6.1f} {m['auc90']:7.3f} {m['duration']:6.1f} "
                  f"{m['attack']:8.1f} {m['aer_share']:9.1f}")
            fh.write(f"{be},{fa},{name},{m['peak']:.3f},{red:.2f},{m['t_peak']:.2f},"
                     f"{m['auc90']:.4f},{m['duration']:.2f},{m['attack']:.2f},"
                     f"{m['aer_share']:.2f}\n")


if __name__ == "__main__":
    main()
