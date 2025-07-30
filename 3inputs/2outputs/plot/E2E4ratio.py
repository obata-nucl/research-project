import matplotlib.pyplot as plt
import numpy as np

from config import BASE_DIR, MODEL_DIR, SAVE_DIR
from spectra import get_expt_energies, get_calc_data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants
# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
MIN_N = 86
MAX_N = 96

SELECTED_LABELS = ["2+_1", "4+_1"]
GRAPH_LABELS = [r"$2^+_1$", r"$4^+_1$"]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def plot_E2E4_ratio(expt_data, calc_data, labels=SELECTED_LABELS):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    markers = ['o', 's']

    calc_colors = ["red", "blue"]
    expt_colors = ["orange", "skyblue"]
    
    # Plot experimental data
    N_array = expt_data["N"].values
    for i in range(len(labels)):
        ax1.plot(N_array, expt_data[labels[i]].values, marker=markers[i], color=expt_colors[i], linestyle='--' ,label=f"Expt. {GRAPH_LABELS[i]}")

    # Plot calculated data
    for i in range(len(labels)):
        ax1.plot(N_array, calc_data[labels[i]].values, marker=markers[i], color=calc_colors[i], label=f"Calc. {GRAPH_LABELS[i]}")
    
    ax1.set_xlabel("Neutron Number")
    ax1.set_ylabel("Energy [MeV]")
    ax1.set_ylim(0,2)
    
    ax2 = ax1.twinx()
    ax2.plot(N_array, expt_data[labels[1]].values / expt_data[labels[0]].values, marker='D', color='grey', linestyle='--', label="Expt. ratio")
    ax2.plot(N_array, calc_data["E(4)/E(2)"].values, marker='D', color="black", label="Calc. ratio")
    ax2.set_ylabel("E(4)/E(2) Ratio")
    ax2.set_ylim(0,4)
    ax2.tick_params(axis='y')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    return fig, ax1, ax2


if __name__ == "__main__":
    expt_data = get_expt_energies()
    calc_data = get_calc_data()

    for model_name, df in calc_data.items():
        save_path = SAVE_DIR / f"{model_name}" / "E2E4ratio.png"
        (SAVE_DIR / model_name).mkdir(parents=True, exist_ok=True)
        fig, ax1, ax2 = plot_E2E4_ratio(expt_data, df, SELECTED_LABELS)
        fig.savefig(save_path)
        plt.close(fig)
