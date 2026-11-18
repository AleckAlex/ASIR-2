"""
paper_figs.py -- regenerate all figures and tables of the revised manuscript.

Uses the model code in asir3.py.  Baseline = winter (poor) ventilation,
zeta = (1.4, 2.9) 1/day; the ventilation measure (zeta x 10) restores
proper ventilation (14, 29) 1/day, the Annex-2 values.

Outputs (in ./paper_out): fig2..fig7 PNG, results.json (all numbers quoted
in the text), table2_peak_reduction.csv, table3_identifiability.csv
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import asir3 as M

OUT = "paper_out"; os.makedirs(OUT, exist_ok=True)
M.ZETA = (1.4, 2.9)                       # winter baseline
ZETA_PROPER = (14.0, 29.0)
FAS = [0.0, 0.2, 0.5, 0.9]
BETAS = [0.3, 0.5]
MEAS = list(M.MITIGATIONS)
COL = {"None": "#55554f", "Barriers": "#2a78d6", "Distancing": "#1baf7a",
       "Ventilation": "#eb6834", "Masking": "#e87ba4", "WFH": "#eda100"}
LS = {"None": "-", "Barriers": "--", "Distancing": ":", "Ventilation": "-.",
      "Masking": (0, (1, 1)), "WFH": "-"}
FA_COL = {0.0: "#b7d3f6", 0.2: "#6da7ec", 0.5: "#2a78d6", 0.9: "#104281"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "savefig.dpi": 300})
R = {}                                     # every number used in the text


def total_I(y):
    return y[M.iI0:M.iI2 + 1].sum(axis=0)


def run(be, fa, name=None, f=None, mult=None, zeta=None):
    """One A-SIR-3 run, calibrated at the winter baseline, seeded on the
    unmitigated growth eigenvector; returns t, y, metrics."""
    d = CAL[(be, fa)]
    y0 = M.initial_state(M.make_params(be * (1 - fa), d))
    if name is not None:
        mb, md, me, mz, budget = M.MITIGATIONS[name]
        mult = (mb, md, me, mz)
        f = M.F_BASE if budget is None else np.array(budget) / 24.0
    f = M.F_BASE if f is None else f
    mult = (1, 1, 1, 1) if mult is None else mult
    old = M.ZETA
    if zeta is not None:
        M.ZETA = zeta
    p = M.make_params(be * (1 - fa), d, f=f, mult=mult)
    M.ZETA = old
    t, y = M.simulate(p, y0)
    return t, y, M.metrics(t, y)


# ---------------------------------------------------------------- calibration
CAL = {}
R["calibration"] = []
for be in BETAS:
    for fa in FAS:
        d = M.calibrate_delta(be, fa)
        CAL[(be, fa)] = d
        p = M.make_params(be * (1 - fa), d)
        R["calibration"].append(dict(beta_eff=be, fa=fa, delta=round(d, 4),
                                     r=round(M.growth_rate(p)[0], 4),
                                     R0=round(M.R0_of(p), 2)))

# ---------------------------------------------------------- Fig 2 validation
def asir1(b, g, d, e, z, i0=0.01, w0=0.01, t_end=100):
    def f(_, y):
        s, i, w = y
        return [-b*s*i - d*s*w, b*s*i + d*s*w - g*i, e*(1-w)*i - z*w]
    t = np.linspace(0, t_end, 4001)
    return t, solve_ivp(f, (0, t_end), [1-i0, i0, w0], t_eval=t, method="Radau",
                        rtol=1e-9, atol=1e-12).y

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
R["validation"] = {}
for k, (b, lab) in enumerate([(0.4, "R₀(SIR) = 2.0"), (0.1, "R₀(SIR) = 0.5")]):
    t, y = asir1(b, 0.2, 0.2, 0.8, 0.05)
    _, ys = asir1(b, 0.2, 0.0, 0.8, 0.05)
    ax[k].plot(t, y[0], color="#2a78d6", lw=1.6, label="S (A-SIR)")
    ax[k].plot(t, y[1], color="#eb6834", lw=1.6, label="I (A-SIR)")
    ax[k].plot(t, y[2], color="#1baf7a", lw=1.6, label="W (A-SIR)")
    ax[k].plot(t, ys[1], color="#eb6834", lw=1.4, ls="--", label="I (SIR)")
    ax[k].set_title(f"β = {b}, γ = 0.2  ({lab})"); ax[k].set_xlabel("Days")
    R["validation"][str(b)] = dict(peak=round(100*y[1].max(), 1),
                                   t_peak=round(t[np.argmax(y[1])], 1),
                                   peak_sir=round(100*ys[1].max(), 1))
ax[0].set_ylabel("Population fraction"); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_validation.png"); plt.close(fig)

# --------------------------------------------------- Fig 3 SIR-3 vs SIR-1
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
R["sir3_vs_sir1"] = {}
for k, be in enumerate(BETAS):
    p = M.make_params(be, 0.0)
    t, y = M.simulate(p, M.initial_state(p))
    t1, I1 = M.sir1(be)
    I3 = total_I(y)
    ax[k].plot(t1, 100*I1, color="#2a78d6", lw=2, label="SIR-1")
    ax[k].plot(t, 100*I3, color="#eb6834", lw=1.6, ls="--", label="SIR-3")
    ax[k].set_xlim(0, 150); ax[k].set_xlabel("Days")
    ax[k].set_title(f"β_eff = {be}, γ = 0.1 (R₀ = {be/0.1:.0f})")
    R["sir3_vs_sir1"][str(be)] = dict(peak3=round(100*I3.max(), 2),
                                      peak1=round(100*I1.max(), 2),
                                      maxdiff_pp=round(100*np.abs(I3-I1).max(), 3))
ax[0].set_ylabel("Infected (%)"); ax[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_sir3_vs_sir1.png"); plt.close(fig)

# --------------------------------------------- mitigation runs (all cases)
RUNS = {}
for be in BETAS:
    for fa in FAS:
        for name in MEAS:
            RUNS[(be, fa, name)] = run(be, fa, name)

def red(be, fa, name):
    return 100 * (1 - RUNS[(be, fa, name)][2]["peak"] / RUNS[(be, fa, "None")][2]["peak"])

R["mitigation"] = []
with open(f"{OUT}/table2_peak_reduction.csv", "w") as fh:
    fh.write("R0,aerosol_pct," + ",".join(MEAS) + "\n")
    for be in BETAS:
        for fa in FAS:
            vals = [RUNS[(be, fa, "None")][2]["peak"]] + [red(be, fa, n) for n in MEAS[1:]]
            fh.write(f"{be/0.1:.0f},{int(fa*100)}," + ",".join(f"{v:.1f}" for v in vals) + "\n")
            for n in MEAS:
                m = RUNS[(be, fa, n)][2]
                R["mitigation"].append(dict(R0=round(be/0.1), fa=fa, measure=n,
                    peak=round(m["peak"], 1), red=round(red(be, fa, n), 1),
                    attack=round(m["attack"], 1), duration=round(m["duration"], 1),
                    aer_share=round(m["aer_share"], 1), t_peak=round(m["t_peak"], 1)))

# Fig 4: grid
fig, ax = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex=True, sharey=True)
for r_, be in enumerate(BETAS):
    for c_, fa in enumerate([0.0, 0.5, 0.9]):
        a = ax[r_, c_]
        for n in MEAS:
            t, y, _ = RUNS[(be, fa, n)]
            a.plot(t, 100*total_I(y), color=COL[n], ls=LS[n], lw=1.5, label=n)
        a.set_xlim(0, 250)
        a.set_title(f"R₀ = {be/0.1:.0f} | {int(fa*100)}% aerosol", fontsize=8.5)
    ax[r_, 0].set_ylabel("Infected (%)")
for a in ax[1]:
    a.set_xlabel("Days")
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=6, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(f"{OUT}/fig4_mitigation_grid.png"); plt.close(fig)

# Fig 5: metrics bars (R0 = 3)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
be = 0.3; x = np.arange(len(MEAS)); w = 0.26
for j, (key, lab) in enumerate([("peak", "Peak infected (%)"),
                                ("attack", "Final attack rate (%)"),
                                ("duration", "Duration above 1% (days)")]):
    for k, fa in enumerate([0.0, 0.5, 0.9]):
        v = [RUNS[(be, fa, n)][2][key] for n in MEAS]
        ax[j].bar(x + (k-1)*w, v, w*0.92, color=FA_COL[fa],
                  label=f"{int(fa*100)}% aerosol")
    ax[j].set_xticks(x); ax[j].set_xticklabels(MEAS, rotation=45, ha="right", fontsize=7.5)
    ax[j].set_title(lab, fontsize=8.5); ax[j].grid(axis="x", visible=False)
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(f"{OUT}/fig5_metrics.png"); plt.close(fig)

# --------------------------------------- Fig 6 aerosol-specific strategies
SUMMER = np.array([4.0, 11.0, 9.0]) / 24.0          # +2 h outdoors
STRAT = {
    "None (winter baseline)":  dict(),
    "Outdoor shift (+2 h)":    dict(f=SUMMER),
    "Home ventilation ζ₁×5":   dict(mult=(1, 1, 1, 1), zeta=(7.0, 2.9)),
    "Work ventilation ζ₂×10":  dict(mult=(1, 1, 1, 1), zeta=(1.4, 29.0)),
    "Home + work ventilation": dict(mult=(1, 1, 1, 1), zeta=(7.0, 29.0)),
    "Outdoor + ventilation":   dict(f=SUMMER, zeta=(7.0, 29.0)),
}
SCOL = ["#55554f", "#1baf7a", "#2a78d6", "#eb6834", "#e87ba4", "#eda100"]
SLS = ["-", "--", ":", "-.", (0, (1, 1)), "-"]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
R["strategies"] = {}
for k, fa in enumerate([0.5, 0.9]):
    base = None
    for j, (lab, kw) in enumerate(STRAT.items()):
        t, y, m = run(0.3, fa, **kw)
        base = m["peak"] if base is None else base
        ax[k].plot(t, 100*total_I(y), color=SCOL[j], ls=SLS[j], lw=1.5, label=lab)
        R["strategies"][f"{fa}|{lab}"] = dict(peak=round(m["peak"], 1),
                                            red=round(100*(1-m["peak"]/base), 1))
    ax[k].set_xlim(0, 250); ax[k].set_xlabel("Days")
    ax[k].set_title(f"R₀ = 3 | {int(fa*100)}% aerosol", fontsize=8.5)
ax[0].set_ylabel("Infected (%)")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.85))
fig.savefig(f"{OUT}/fig6_strategies.png"); plt.close(fig)

# ------------------------------------- Fig 7 identifiability gap (R0 = 3)
rng = np.random.default_rng(2026)
SIG, WIN = 0.002, 60
SIR_MULT = {"Barriers": 0.4, "Distancing": 0.5, "Ventilation": 1.0,
            "Masking": 0.1, "WFH": 0.2}

def sir_curve(b, g, i0, t):
    def f(_, y):
        s, i = y
        return [-b*s*i, b*s*i - g*i]
    return solve_ivp(f, (0, t[-1]), [1-i0, i0], t_eval=t, method="LSODA",
                     rtol=1e-9, atol=1e-12).y[1]

R["identifiability"] = []
fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.1])
gap_rows = []
for k, fa in enumerate(FAS):
    t, y, _ = RUNS[(0.3, fa, "None")]
    tf = np.arange(0, WIN + 1, 1.0)
    data = np.interp(tf, t, total_I(y)) + rng.normal(0, SIG, tf.size)
    fits = [least_squares(lambda p: sir_curve(p[0], p[1], M.I_INIT, tf) - data, x0,
                          bounds=([0.01, 0.01], [3, 3]))
            for x0 in ([0.3, 0.1], [0.5, 0.2], [0.2, 0.05])]
    fit = min(fits, key=lambda r: r.cost)
    b_hat, g_hat = fit.x
    rms = np.sqrt(2*fit.cost/tf.size)
    tt = np.linspace(0, M.T_END, 1441)
    base_sir = sir_curve(b_hat, g_hat, M.I_INIT, tt).max()
    a = fig.add_subplot(gs[0, k])
    a.plot(tf, 100*data, ".", ms=2.5, color="#898781", label="data (A-SIR-3)")
    a.plot(tt, 100*np.interp(tt, t, total_I(y)), color="#eb6834", lw=1.5, label="truth")
    a.plot(tt, 100*sir_curve(b_hat, g_hat, M.I_INIT, tt), color="#2a78d6", lw=1.3,
           ls="--", label="SIR-1 fit")
    a.set_xlim(0, 120); a.set_xlabel("Days", fontsize=8)
    a.set_title(f"{int(fa*100)}% aerosol\nR̂₀ = {b_hat/g_hat:.2f}", fontsize=8)
    if k == 0:
        a.set_ylabel("Infected (%)"); a.legend(frameon=False, fontsize=6.5)
    row = dict(fa=fa, beta_hat=round(b_hat, 4), gamma_hat=round(g_hat, 4),
               R0_hat=round(b_hat/g_hat, 3), rms_resid_pct=round(100*rms, 3))
    for n, mb in SIR_MULT.items():
        pk = sir_curve(b_hat*mb, g_hat, M.I_INIT, tt).max()
        pred = 100*(1 - pk/base_sir)
        true = red(0.3, fa, n)
        row[n] = dict(sir=round(pred, 1), true=round(true, 1), gap=round(pred-true, 1))
    R["identifiability"].append(row)
a = fig.add_subplot(gs[1, :])
xx = np.arange(len(SIR_MULT)); w = 0.2
for k, fa in enumerate(FAS):
    gaps = [R["identifiability"][k][n]["gap"] for n in SIR_MULT]
    bars = a.bar(xx + (k-1.5)*w, gaps, w*0.92, color=FA_COL[fa], label=f"{int(fa*100)}% aerosol")
a.axhline(0, color="#55554f", lw=0.8)
a.set_xticks(xx); a.set_xticklabels(list(SIR_MULT))
a.set_ylabel("SIR-predicted minus true\npeak reduction (pp)")
a.legend(frameon=False, ncol=4, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, 1.0)); a.grid(axis="x", visible=False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig7_identifiability.png"); plt.close(fig)

with open(f"{OUT}/table3_identifiability.csv", "w") as fh:
    fh.write("aerosol_pct,R0_hat,measure,SIR_pred,true,gap\n")
    for row in R["identifiability"]:
        for n in SIR_MULT:
            fh.write(f"{int(row['fa']*100)},{row['R0_hat']},{n},{row[n]['sir']},"
                     f"{row[n]['true']},{row[n]['gap']}\n")

# proper-ventilation sensitivity (baseline zeta = proper; ventilation x10 beyond)
R["proper_sensitivity"] = {}
oldZ = M.ZETA; M.ZETA = ZETA_PROPER
for be in BETAS:
    for fa in [0.5, 0.9]:
        d = M.calibrate_delta(be, fa)
        y0 = M.initial_state(M.make_params(be*(1-fa), d))
        out = {}
        for n in MEAS:
            mb, md, me, mz, bud = M.MITIGATIONS[n]
            f = M.F_BASE if bud is None else np.array(bud)/24
            t, y = M.simulate(M.make_params(be*(1-fa), d, f=f, mult=(mb, md, me, mz)), y0)
            out[n] = M.metrics(t, y)["peak"]
        R["proper_sensitivity"][f"{be}|{fa}"] = {n: round(100*(1-out[n]/out["None"]), 1)
                                                for n in MEAS[1:]}
M.ZETA = oldZ

json.dump(R, open(f"{OUT}/results.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps({k: R[k] for k in ("calibration", "validation", "sir3_vs_sir1",
                                    "strategies", "identifiability", "proper_sensitivity")},
                 indent=1, ensure_ascii=False))
