import matplotlib.pyplot as plt
import numpy as np

from config import BASE_DIR, MODEL_DIR, SAVE_DIR
from spectra import get_calc_data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants

# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
MIN_N = 86
MAX_N = 96
N_NU = np.array([2, 3, 4, 5, 6, 7])
BETA_F = np.arange(-0.5, 0.5, 0.1)

PARAMS_LABEL = ["eps", "kappa"]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def IBM2_energy(output, n_nu, beta_f):
    # IBM parameters
    eps = output[0]
    kappa = output[1]

    beta_b = beta_f * 5
    beta2 = beta_b * beta_b
    N = 6 + n_nu
    deno = 1 / (1 + beta2)

    # Calculate the IBM2 energy
    IBM2_energy = (
        eps * N * beta2 * deno
    ) + (
        6 * n_nu * kappa * beta2 * deno * deno * (4 + 1.3 * 2 * np.sqrt(2 / 7) * beta_b + 0.8 * beta2 / 7)
    )

    return IBM2_energy



def plot_PES(N:int, n_nu: int, params: np.ndarray, beta_f: np.ndarray, gamma_f: np.ndarray = None):
    if gamma_f is None:
        Energy = IBM2_energy(params, n_nu, beta_f)
        Energy -= Energy.min()

        # Find the global minimum
        min_idx_global = np.argmin(Energy)
        min_beta_val = beta_f[min_idx_global]

        # Check if the minimum is very close to zero
        if np.isclose(min_beta_val, 0):
            # If so, find the index corresponding to beta_f = 0
            zero_idx = np.argmin(np.abs(beta_f))
            min_idx_to_plot = zero_idx
        else:
            # Otherwise, use the actual minimum for beta_f >= 0
            mask = beta_f >= 0
            energy_pos = Energy[mask]
            if len(energy_pos) > 0:
                min_idx_local = np.argmin(energy_pos)
                min_idx_to_plot = np.where(mask)[0][min_idx_local]
            else: # If no beta_f >= 0, just use global minimum
                min_idx_to_plot = min_idx_global


        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(beta_f, Energy, linestyle='-', color="black")
        ax.plot(beta_f[min_idx_to_plot], Energy[min_idx_to_plot], 'ro', markersize=8)
        ax.set_xlabel(r"$\beta_f$")
        ax.set_ylabel("Energy [MeV]")
        ax.set_title(f"{N+62}Sm")

        return fig, ax



if __name__ == "__main__":
    for model_name, df in get_calc_data().items():
        save_dir = SAVE_DIR / model_name / "PES"
        save_dir.mkdir(parents=True, exist_ok=True)
        inputs = df[PARAMS_LABEL].to_numpy()
        N_array = np.arange(MIN_N, MAX_N + 1, 2)
        for i in range(len(N_array)):
            fig, ax = plot_PES(N_array[i], N_NU[i], inputs[i], BETA_F)
            fig.savefig(save_dir / f"{N_array[i]}.png")
            plt.close(fig)