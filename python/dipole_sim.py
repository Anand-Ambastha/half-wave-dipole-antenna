#!/usr/bin/env python3
"""
Half-wave dipole: component-level Python simulation (method of moments)
=========================================================================
Companion to the CST model (dipole_2p4GHz.bas).  Components:

  1. Geometry        : Dipole dataclass (length, radius, gap, material)
  2. MoM solver      : thin-wire mixed-potential EFIE, triangle basis, Galerkin
  3. Feed / port     : delta-gap source -> Zin, reflection coefficient S11, VSWR
  4. Loss model      : surface-resistance (skin effect) ohmic loss -> efficiency
  5. Far field       : radiation pattern, directivity, gain, HPBW
  6. Analysis        : frequency sweep, length sweep, radius sweep, optimiser
  7. Plots / report  : PNG figures + results.csv + summary.txt

Usage:   python3 dipole_sim.py            (full run, writes ./dipole_results/)
         python3 dipole_sim.py --converge (mesh-convergence check only)
Needs:   numpy, scipy, matplotlib
"""
import argparse
import os
from dataclasses import dataclass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

C0 = 299792458.0
MU0 = 4e-7 * np.pi
EPS0 = 1.0 / (MU0 * C0**2)
ETA0 = MU0 * C0
Z0 = 50.0


# ---------------------------------------------------------------- 1. geometry
@dataclass
class Dipole:
    L: float = 58.7e-3          # tip-to-tip length [m]
    a: float = 1.0e-3           # wire radius [m]
    sigma: float = np.inf       # conductivity [S/m]  (inf = PEC; copper = 5.8e7)
    gap: float = 2.0e-3         # feed gap width [m] (same as the CST discrete port)

    def with_(self, **kw):
        d = dict(L=self.L, a=self.a, sigma=self.sigma, gap=self.gap)
        d.update(kw)
        return Dipole(**d)


# ---------------------------------------------------------------- 2. MoM solver
def _segments_per_wire(dp, n_min=20):
    """Even number of segments, segment length ~ max(gap/2, radius)."""
    h = max(dp.gap / 2.0, dp.a)
    n = max(int(round(dp.L / h)), n_min)
    return n + (n % 2)           # even so a node sits at the feed point


def solve_current(dp, f, n_seg=None, nq=6):
    """Solve for the current on the dipole with a 1 V delta-gap at the centre.

    Mixed-potential EFIE with reduced (thin-wire) kernel
        G(z,z') = exp(-jkR)/(4 pi R),  R = sqrt((z-z')^2 + a^2)
    and triangle basis / testing functions (Galerkin).  The current is zero at
    both wire ends by construction (free-end boundary condition).

    Returns (z_nodes, I_nodes, Zin)
    """
    k = 2 * np.pi * f / C0
    w = 2 * np.pi * f
    N = n_seg if n_seg is not None else _segments_per_wire(dp)
    h = dp.L / N
    nodes = np.linspace(-dp.L / 2, dp.L / 2, N + 1)
    M = N - 1                                   # interior nodes = unknowns

    # Gauss-Legendre points on one segment
    xg, wg = np.polynomial.legendre.leggauss(nq)
    xs = (xg + 1) / 2 * h                       # offsets in [0,h]
    ws = wg * h / 2
    # triangle values at the quadrature points: left half rising, right falling
    rise = xs / h
    fall = 1 - rise

    # Pre-compute segment-to-segment interaction of the 4 products
    # (basis n on segment s', testing m on segment s); then assemble.
    seg_start = nodes[:-1]
    Zmat = np.zeros((M, M), dtype=complex)

    # Segment-pair kernel integrals: I_pq[s,s'] for p,q in {rise,fall}
    # kernel matrix between all quadrature points of segs s and s'
    zq = seg_start[:, None] + xs[None, :]       # (N, nq)
    # G between every pair of quadrature points
    dz = zq[:, None, :, None] - zq[None, :, None, :]       # (N,N,nq,nq)
    R = np.sqrt(dz**2 + dp.a**2)
    G = np.exp(-1j * k * R) / (4 * np.pi * R)
    wq = ws[None, None, :, None] * ws[None, None, None, :]

    # Kernel integrals between segment s (testing) and segment t (basis) for each
    # combination of rising/falling triangle halves. The quadrature weights are
    # already contained in wq.
    A_rr = np.einsum("stij,i,j->st", G * wq, rise, rise)   # test rise, basis rise
    A_rf = np.einsum("stij,i,j->st", G * wq, rise, fall)
    A_fr = np.einsum("stij,i,j->st", G * wq, fall, rise)
    A_ff = np.einsum("stij,i,j->st", G * wq, fall, fall)
    S = np.einsum("stij->st", G * wq)                       # derivative products

    # Triangle at node n (1..N-1): rising on segment n-1, falling on segment n.
    # The testing functions use the same triangles (Galerkin).
    for m in range(1, N):
        for n in range(1, N):
            sm_l, sm_r = m - 1, m          # segments of testing function m
            sn_l, sn_r = n - 1, n          # segments of basis function n
            vec = (A_rr[sm_l, sn_l] + A_rf[sm_l, sn_r]
                   + A_fr[sm_r, sn_l] + A_ff[sm_r, sn_r])
            # derivative of triangle: +1/h on left segment, -1/h on right
            sc = (S[sm_l, sn_l] * (+1) * (+1) + S[sm_l, sn_r] * (+1) * (-1)
                  + S[sm_r, sn_l] * (-1) * (+1) + S[sm_r, sn_r] * (-1) * (-1)) / h**2
            Zmat[m - 1, n - 1] = 1j * w * MU0 * vec + sc / (1j * w * EPS0)

    # Feed: 1 V applied as a uniform field V/gap across the feed gap, exactly
    # what a CST discrete port does.  V_m = (V/gap) * integral of T_m over gap.
    v = np.zeros(M, dtype=complex)
    feed = N // 2 - 1                      # index of centre node in unknowns
    zg = np.linspace(-dp.gap / 2, dp.gap / 2, 401)
    for m in range(1, N):
        tm = np.clip(1 - np.abs(zg - nodes[m]) / h, 0, None)
        v[m - 1] = np.trapezoid(tm, zg) / dp.gap
    I = np.linalg.solve(Zmat, v)
    Zin = 1.0 / I[feed]

    I_nodes = np.concatenate([[0], I, [0]])
    return nodes, I_nodes, Zin


# ---------------------------------------------------------------- 3. feed/port
def s11_from_z(zin, z0=Z0):
    return (zin - z0) / (zin + z0)


def vswr_from_s11(s11):
    g = np.abs(s11)
    return (1 + g) / (1 - g)


def freq_sweep(dp, freqs, n_seg=None):
    z = np.array([solve_current(dp, f, n_seg)[2] for f in freqs])
    s11 = s11_from_z(z)
    return z, s11


# ---------------------------------------------------------------- 4. loss model
def ohmic_loss_resistance(dp, f, z, I):
    """Equivalent series resistance of conductor loss referred to the feed."""
    if np.isinf(dp.sigma):
        return 0.0
    Rs = np.sqrt(np.pi * f * MU0 / dp.sigma)           # surface resistance
    Rl = Rs / (2 * np.pi * dp.a)                       # ohm/m
    # P_loss = 0.5 * Rl * integral |I|^2 dz ; refer to feed current
    p = 0.5 * Rl * np.trapezoid(np.abs(I) ** 2, z)
    i_feed = I[len(I) // 2]
    return 2 * p / np.abs(i_feed) ** 2


# ---------------------------------------------------------------- 5. far field
def farfield(dp, f, z, I, theta):
    """E_theta pattern (unnormalised) from the current distribution."""
    k = 2 * np.pi * f / C0
    th = np.asarray(theta)
    # integral of I(z') exp(j k z' cos(theta)) dz' with piecewise-linear current
    integrand = I[None, :] * np.exp(1j * k * np.outer(np.cos(th), z))
    F = np.trapezoid(integrand, z, axis=1)
    return np.sin(th) * F


def directivity_and_hpbw(dp, f, z, I, ntheta=2001):
    th = np.linspace(1e-6, np.pi - 1e-6, ntheta)
    E = farfield(dp, f, z, I, th)
    U = np.abs(E) ** 2
    prad = 2 * np.pi * np.trapezoid(U * np.sin(th), th)
    D = 4 * np.pi * U / prad                           # directivity vs theta
    Dmax = D.max()
    # E-plane HPBW (around theta = 90 deg)
    half = Dmax / 2
    i_max = np.argmax(D)
    left = np.where(D[:i_max] < half)[0]
    right = np.where(D[i_max:] < half)[0]
    th1 = th[left[-1]] if left.size else th[0]
    th2 = th[i_max + right[0]] if right.size else th[-1]
    hpbw = np.degrees(th2 - th1)
    return th, D, Dmax, hpbw, prad


def radiation_efficiency(dp, f, z, I, zin):
    Rl = ohmic_loss_resistance(dp, f, z, I)
    Rin = zin.real
    return (Rin - Rl) / Rin if Rin > 0 else np.nan


# ---------------------------------------------------------------- 6. analysis
def resonant_length(f, a, lo=0.40, hi=0.50, sigma=np.inf):
    """Length (m) where Im(Zin)=0 at frequency f for radius a."""
    lam = C0 / f

    def x(Lr):
        return solve_current(Dipole(L=Lr * lam, a=a, sigma=sigma), f)[2].imag

    # scan for the first sign change of the reactance, then refine with brentq
    grid = np.linspace(0.30, 0.52, 23)
    xs = [x(g) for g in grid]
    for i in range(len(grid) - 1):
        if xs[i] < 0 <= xs[i + 1]:
            return brentq(x, grid[i], grid[i + 1], xtol=1e-5) * lam
    raise RuntimeError("no resonance found in 0.30-0.52 lambda")


def bandwidth(freqs, s11, thresh_db=-10.0):
    """Contiguous -10 dB band around the S11 minimum."""
    s = 20 * np.log10(np.abs(s11))
    i0 = np.argmin(s)
    if s[i0] > thresh_db:
        return None, None, None
    lo = i0
    while lo > 0 and s[lo - 1] <= thresh_db:
        lo -= 1
    hi = i0
    while hi < len(s) - 1 and s[hi + 1] <= thresh_db:
        hi += 1
    return freqs[lo], freqs[hi], freqs[i0]


# ---------------------------------------------------------------- 7. plotting
def smith_grid(ax):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), "k", lw=0.8)
    for r in [0.2, 0.5, 1, 2, 5]:
        c = r / (1 + r)
        rad = 1 / (1 + r)
        ax.plot(c + rad * np.cos(t), rad * np.sin(t), color="0.75", lw=0.6)
    for x in [0.2, 0.5, 1, 2, 5]:
        for sgn in (1, -1):
            xx = sgn * x
            cc = 1j * 1 / xx + 1
            rad = 1 / abs(xx)
            tt = np.linspace(0, 2 * np.pi, 800)
            pts = cc + rad * np.exp(1j * tt)
            m = np.abs(pts) <= 1.0001
            ax.plot(pts.real[m], pts.imag[m], color="0.75", lw=0.6)
    ax.axhline(0, color="0.75", lw=0.6)
    ax.set_aspect("equal")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.axis("off")


def make_plots(out, dp, f0, freqs, z, s11, zc, Ic, th, D, lens, zl, rads, zr, res_len):
    lam = C0 / f0

    # current distribution
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(zc * 1e3, np.abs(Ic) * 1e3, label="|I|")
    ax.plot(zc * 1e3, Ic.real * 1e3, "--", label="Re(I)")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("current for 1 V feed [mA]")
    ax.set_title(f"Current distribution @ {f0/1e9:.2f} GHz (L={dp.L*1e3:.1f} mm)")
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "current.png"), dpi=150)
    plt.close(fig)

    # S11, VSWR
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(freqs / 1e9, 20 * np.log10(np.abs(s11)))
    ax[0].axhline(-10, color="r", ls="--", lw=0.8)
    ax[0].set_xlabel("Frequency [GHz]")
    ax[0].set_ylabel("S11 [dB]")
    ax[0].set_title("Return loss (50 ohm)")
    ax[0].grid(alpha=.3)
    ax[1].plot(freqs / 1e9, vswr_from_s11(s11))
    ax[1].axhline(2, color="r", ls="--", lw=0.8)
    ax[1].set_ylim(1, 10)
    ax[1].set_xlabel("Frequency [GHz]")
    ax[1].set_ylabel("VSWR")
    ax[1].set_title("VSWR")
    ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "s11_vswr.png"), dpi=150)
    plt.close(fig)

    # input impedance
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(freqs / 1e9, z.real, label="R")
    ax.plot(freqs / 1e9, z.imag, label="X")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("Zin [ohm]")
    ax.set_title("Input impedance")
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "zin.png"), dpi=150)
    plt.close(fig)

    # Smith chart
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    smith_grid(ax)
    ax.plot(s11.real, s11.imag, "b", lw=1.5)
    ax.plot(s11.real[0], s11.imag[0], "go")
    ax.plot(s11.real[-1], s11.imag[-1], "rs")
    ax.set_title("Smith chart (green start, red end)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "smith.png"), dpi=150)
    plt.close(fig)

    # radiation patterns (E-plane polar, H-plane polar)
    Dn = D / D.max()
    Ddb = 10 * np.log10(np.maximum(D, 1e-6))
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(121, projection="polar")
    phi = np.concatenate([th, 2 * np.pi - th[::-1]])
    r = np.concatenate([Ddb, Ddb[::-1]])
    ax1.plot(phi, np.maximum(r, -20), lw=1.5)
    ax1.set_theta_zero_location("N")
    ax1.set_theta_direction(-1)
    ax1.set_rlim(-20, Ddb.max() + 1)
    ax1.set_title("E-plane (vertical, dBi)")
    ax2 = fig.add_subplot(122, projection="polar")
    ph = np.linspace(0, 2 * np.pi, 361)
    ax2.plot(ph, np.full_like(ph, Ddb.max()), lw=1.5)
    ax2.set_rlim(-20, Ddb.max() + 1)
    ax2.set_title("H-plane (azimuth, dBi), omnidirectional")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pattern_polar.png"), dpi=150)
    plt.close(fig)

    # 3D pattern
    ph = np.linspace(0, 2 * np.pi, 91)
    TH, PH = np.meshgrid(th, ph)
    Rr = np.tile(Dn, (len(ph), 1))
    X = Rr * np.sin(TH) * np.cos(PH)
    Y = Rr * np.sin(TH) * np.sin(PH)
    Zz = Rr * np.cos(TH)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Zz, facecolors=plt.cm.viridis(Rr), rstride=2, cstride=2,
                    linewidth=0, antialiased=False)
    ax.set_title("3D pattern (normalised)")
    ax.set_box_aspect((1, 1, 1))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pattern_3d.png"), dpi=150)
    plt.close(fig)

    # length sweep
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lens / lam, zl.real, label="R")
    ax.plot(lens / lam, zl.imag, label="X")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(res_len / lam, color="g", ls="--", label=f"resonance {res_len*1e3:.1f} mm")
    ax.set_xlabel("L / lambda")
    ax.set_ylabel("Zin [ohm]")
    ax.set_title(f"Length sweep @ {f0/1e9:.1f} GHz, a = {dp.a*1e3:.1f} mm")
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "length_sweep.png"), dpi=150)
    plt.close(fig)

    # radius sweep S11
    fig, ax = plt.subplots(figsize=(6, 4))
    for a, (fr, s) in rads.items():
        ax.plot(fr / 1e9, 20 * np.log10(np.abs(s)), label=f"a = {a*1e3:.2f} mm")
    ax.axhline(-10, color="r", ls="--", lw=0.8)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("S11 [dB]")
    ax.set_title("Radius sweep (each length re-tuned to 2.4 GHz)")
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "radius_sweep.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- convergence
def convergence_check(f=2.4e9):
    lam = C0 / f
    dp = Dipole(L=0.4 * lam + 0 * lam, a=1e-3)   # placeholder replaced below
    print("Mesh convergence of the gap-fed model (a = 1 mm, gap = 2 mm):")
    for Lr in (0.47, 0.50):
        dp = Dipole(L=Lr * lam, a=1e-3)
        print(f"  L = {Lr} lambda = {dp.L*1e3:.1f} mm")
        for h in (2.0e-3, 1.5e-3, 1.0e-3, 0.75e-3):
            N = int(round(dp.L / h)); N += N % 2
            _, _, zin = solve_current(dp, f, n_seg=N)
            print(f"    N={N:3d} h={dp.L/N*1e3:4.2f} mm  Zin = {zin.real:6.2f} {zin.imag:+7.2f}j")
    print("Reference: thin lambda/2 dipole ~ 73 + j42 ohm (a -> 0)")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--converge", action="store_true")
    ap.add_argument("--f0", type=float, default=2.4e9)
    ap.add_argument("--a", type=float, default=1.0e-3)
    ap.add_argument("--out", default="dipole_results")
    args = ap.parse_args()

    if args.converge:
        convergence_check(args.f0)
        return

    os.makedirs(args.out, exist_ok=True)
    f0 = args.f0
    lam = C0 / f0
    a = args.a

    print("1/6 Finding resonant length ...")
    L_res = resonant_length(f0, a)
    print(f"    resonant length = {L_res*1e3:.2f} mm = {L_res/lam:.3f} lambda")
    dp = Dipole(L=L_res, a=a)

    print("2/6 Solving current at f0 ...")
    zc, Ic, zin0 = solve_current(dp, f0)
    print(f"    Zin(f0) = {zin0.real:.1f} {zin0.imag:+.1f}j ohm")

    print("3/6 Frequency sweep ...")
    freqs = np.linspace(1.5e9, 3.5e9, 81)
    z, s11 = freq_sweep(dp, freqs)
    flo, fhi, fmin = bandwidth(freqs, s11)

    print("4/6 Far field ...")
    th, D, Dmax, hpbw, prad = directivity_and_hpbw(dp, f0, zc, Ic)
    eff_pec = 1.0
    cu = dp.with_(sigma=5.8e7)
    zc2, Ic2, zin_cu = solve_current(cu, f0)
    eff_cu = radiation_efficiency(cu, f0, zc2, Ic2, zin_cu)

    print("5/6 Length and radius sweeps ...")
    lens = np.linspace(0.40, 0.50, 21) * lam
    zl = np.array([solve_current(dp.with_(L=Lx), f0)[2] for Lx in lens])
    rads = {}
    rad_rows = []
    for ar in [0.25e-3, 0.5e-3, 1.0e-3]:
        Lr = resonant_length(f0, ar)
        fr = np.linspace(1.8e9, 3.0e9, 49)
        zr, sr = freq_sweep(Dipole(L=Lr, a=ar), fr)
        rads[ar] = (fr, sr)
        b = bandwidth(fr, sr)
        bw = (b[1] - b[0]) / f0 * 100 if b[0] else 0.0
        rad_rows.append((ar, Lr, bw))

    print("6/6 Plots and report ...")
    make_plots(args.out, dp, f0, freqs, z, s11, zc, Ic, th, D, lens, zl, rads, None, L_res)

    # CSV of the sweep
    with open(os.path.join(args.out, "results.csv"), "w") as fh:
        fh.write("f_GHz,R_ohm,X_ohm,S11_dB,VSWR\n")
        for fq, zz, ss in zip(freqs, z, s11):
            fh.write(f"{fq/1e9:.4f},{zz.real:.3f},{zz.imag:.3f},"
                     f"{20*np.log10(abs(ss)):.3f},{vswr_from_s11(ss):.3f}\n")

    s11_f0 = s11_from_z(zin0)
    lines = [
        "HALF-WAVE DIPOLE: PYTHON (MoM) RESULTS",
        "=" * 44,
        f"Design frequency        : {f0/1e9:.3f} GHz  (lambda = {lam*1e3:.2f} mm)",
        f"Wire radius             : {a*1e3:.2f} mm  (a/lambda = {a/lam:.4f})",
        f"Resonant length (X = 0) : {L_res*1e3:.2f} mm  = {L_res/lam:.3f} lambda",
        f"Zin at f0               : {zin0.real:.1f} {zin0.imag:+.1f}j ohm",
        f"S11 at f0 (50 ohm)      : {20*np.log10(abs(s11_f0)):.1f} dB",
        f"VSWR at f0              : {vswr_from_s11(s11_f0):.2f}",
        f"Max directivity         : {10*np.log10(Dmax):.2f} dBi",
        f"E-plane HPBW            : {hpbw:.1f} deg",
        f"Efficiency (PEC)        : {eff_pec*100:.1f} %",
        f"Efficiency (copper)     : {eff_cu*100:.2f} %  -> gain {10*np.log10(Dmax*eff_cu):.2f} dBi",
    ]
    if flo:
        lines.append(f"-10 dB band             : {flo/1e9:.3f}-{fhi/1e9:.3f} GHz "
                     f"({(fhi-flo)/f0*100:.1f} % ), min S11 at {fmin/1e9:.3f} GHz")
    else:
        lines.append("-10 dB band             : not reached against 50 ohm")
    lines.append("")
    lines.append("Radius study (length re-tuned to resonance at f0):")
    for ar, Lr, bw in rad_rows:
        lines.append(f"  a = {ar*1e3:.2f} mm -> L = {Lr*1e3:.2f} mm, -10 dB BW = {bw:.1f} %")
    lines.append("")
    lines.append("Figures: current, s11_vswr, zin, smith, pattern_polar, pattern_3d,")
    lines.append("         length_sweep, radius_sweep  (PNG, in the output folder)")
    text = "\n".join(lines)
    print("\n" + text)
    with open(os.path.join(args.out, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
