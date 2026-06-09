import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

# Constants & Params
m2, rh, eps, rmax, n_eval = -2.0, 1.0, 1e-6, 30.0, 4000
rh3 = rh**3
T_h = 3.0 * rh / (4.0 * np.pi)
t_ev = np.linspace(rh + eps, rmax, n_eval)

# --- PRECOMPUTE ASYMPTOTIC MATRICES FOR MASSIVE SPEEDUP ---
# Because t_ev is fixed, we can compute the pseudo-inverse of the fitting matrices ONCE.
mask = t_ev > 0.7 * rmax
rt = t_ev[mask]
inv_rt = 1.0 / rt
A_phi_pinv = np.linalg.pinv(np.column_stack([np.ones_like(rt), inv_rt]))
A_psi_pinv = np.linalg.pinv(np.column_stack([inv_rt, inv_rt**2]))

def odes(r, y):
    """System of 4 first-order ODEs for psi, psi', phi, phi'."""
    psi, psip, phi, phip = y
    inv_r = 1.0 / r
    f = r**2 - rh3 * inv_r
    inv_f = 1.0 / f
    fp = 2.0 * r + rh3 * inv_r**2

    return [
        psip,
        -(fp * inv_f + 2.0 * inv_r) * psip - (phi**2 * inv_f**2 - m2 * inv_f) * psi,
        phip,
        -2.0 * inv_r * phip + 2.0 * psi**2 * phi * inv_f,
    ]

def get_asymp(sol):
    """Extract boundary coefficients using fast matrix multiplication."""
    y = sol.y[:, mask]
    mu, n_rho = A_phi_pinv @ y[2]
    psi1, psi2 = A_psi_pinv @ y[0]
    return mu, -n_rho, psi1, psi2

def solve_background(psi_h, E_h):
    """Solve the IVP from the horizon to rmax."""
    y0 = [
        psi_h + (m2 * psi_h / (3.0 * rh)) * eps,
        m2 * psi_h / (3.0 * rh),
        E_h * eps,
        E_h,
    ]
    return solve_ivp(
        odes, (rh + eps, rmax), y0,
        t_eval=t_ev, method="Radau", rtol=1e-7, atol=1e-9
    )

def main():
    results = []
    E_scan = np.linspace(0.01, 30.0, 100)

    print("Running simulations to find source-free solutions...")

    for psi_h in np.linspace(0.05, 2.5, 36):
        source_values = []
        for E in E_scan:
            sol_test = solve_background(psi_h, E)
            source_values.append(get_asymp(sol_test)[2])  # Index 2 is psi1

            if len(source_values) > 1 and source_values[-2] * source_values[-1] < 0:
                root = root_scalar(
                    lambda x: get_asymp(solve_background(psi_h, x))[2],
                    bracket=[E_scan[len(source_values) - 2], E], method="brentq"
                )
                sol = solve_background(psi_h, root.root)
                mu, rho, psi1, psi2 = get_asymp(sol)

                if rho > 0:
                    results.append({
                        "T_hat": T_h / np.sqrt(rho),
                        "O2_hat": (np.sqrt(2.0) * abs(psi2)) / rho,
                        "sol": sol, "mu": mu, "rho": rho, "psi1": psi1, "psi2": psi2
                    })
                break

    if not results:
        print("No valid results found. Adjust scan ranges or solver parameters.")
        return

    results.sort(key=lambda x: x["T_hat"])
    T_hat = np.array([r["T_hat"] for r in results])
    O2_hat = np.array([r["O2_hat"] for r in results])

    # Estimate Tc using the scaled condensate
    mask_Tc = (O2_hat > 0.02 * O2_hat.max()) & (O2_hat < 0.35 * O2_hat.max())
    if mask_Tc.sum() < 2:
        raise RuntimeError("Not enough points near Tc to fit Tc. Adjust scan range or mask.")

    p = np.polyfit(T_hat[mask_Tc], O2_hat[mask_Tc] ** 2, 1)
    Tc_hat = -p[1] / p[0]

    print(f"Estimated Tc/sqrt(rho) = {Tc_hat:.6f}\nReference value is approximately 0.118 for this model.")

    # --- 1. Profiles ---
    profile_indices = sorted(set([0, len(results) // 3, 2 * len(results) // 3, len(results) - 1]))
    for idx, title, ylab in [(0, "Scalar Field", r"$\psi(r)$"), (2, "Gauge Field", r"$\phi(r)$")]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for i in profile_indices:
            ax.plot(results[i]["sol"].t, results[i]["sol"].y[idx], label=rf"$T/T_c={results[i]['T_hat'] / Tc_hat:.3f}$")
        ax.set(title=title, xlabel=r"$r$", ylabel=ylab)
        ax.legend(); ax.grid(True, ls="--")
        plt.tight_layout(); plt.show()

    # --- 2. Phase transition ---
    plt.figure(figsize=(8, 5))
    plt.plot(T_hat / Tc_hat, np.sqrt(O2_hat / Tc_hat ** 2), "o-")
    plt.gca().set(xlabel=r"$T/T_c$", ylabel=r"$\sqrt{\langle \mathcal{O}_2 \rangle}/T_c$", title="Phase Transition")
    plt.grid(True, ls="--")
    plt.tight_layout(); plt.show()

    # --- 3. Asymptotic checks ---
    s_sol = results[len(results) // 2]["sol"]
    mu, rho_s, psi1_s, psi2_s = get_asymp(s_sol)
    r, psi, phi = s_sol.t, s_sol.y[0], s_sol.y[2]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(r, r * psi, lw=2)
    ax1.axhline(0, color="k", ls="--", alpha=0.6)
    ax1.set(xlabel=r"$r$", ylabel=r"$r\psi(r)$", title=rf"Source check: $\psi_1 \approx {psi1_s:.2e}$")
    ax1.grid(True, ls="--")

    ax2.plot(r, r**2 * psi, lw=2)
    ax2.axhline(psi2_s, color="k", ls="--", label=rf"$\psi_2 \approx {psi2_s:.3f}$")
    ax2.set(xlabel=r"$r$", ylabel=r"$r^2\psi(r)$", title="Normalisable Mode")
    ax2.legend(); ax2.grid(True, ls="--")
    plt.tight_layout(); plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(1.0 / r, phi, lw=2, label=r"Solution $\phi(r)$")
    plt.plot(1.0 / r, mu - rho_s / r, "k--", label=r"Fit $\mu - \rho/r$", alpha=0.8)
    plt.plot(1.0 / rh, 0, "ro", label="Horizon", ms=8)
    plt.gca().invert_xaxis()
    plt.gca().set(xlabel=r"$1/r$", ylabel=r"$\phi(r)$", title="Gauge Asymptotics")
    plt.legend(); plt.grid(True, ls="--")
    plt.tight_layout(); plt.show()

if __name__ == "__main__":
    main()