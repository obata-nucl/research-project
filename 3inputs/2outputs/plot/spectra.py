import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import BASE_DIR, MODEL_DIR, SAVE_DIR

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants
# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
MIN_N = 86
MAX_N = 96

SELECTED_LABELS = ["2+_1", "4+_1", "6+_1", "0+_2"]
GRAPH_LABELS = [r"$2^+_1$", r"$4^+_1$", r"$6^+_1$", r"$0^+_2$"]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def get_expt_energies() -> pd.DataFrame:
    expt_path = BASE_DIR / "data" / "62" / "expt.csv"
    df = pd.read_csv(expt_path)
    filtered_df = df[(df["N"] >= MIN_N) & (df["N"] <= MAX_N)]
    return filtered_df

def get_calc_data() -> dict[str, pd.DataFrame]:
    data_list = {}
    for each_file in (MODEL_DIR / "results_data" / "62").iterdir():
        if each_file.suffix == ".csv":
            df = pd.read_csv(each_file)
            data_list[each_file.stem] = df
    return data_list



def plot_energies(N_array, energies, labels=GRAPH_LABELS):
    fig, ax = plt.subplots(figsize=(10, 6))
    markers = ['o', 's', '^', 'D']
    for i in range(len(energies)):
        ax.plot(N_array, energies[i], marker=markers[i], label=labels[i])
    ax.set_xlabel("Neutron Number")
    ax.set_ylim(0, 2)
    ax.legend()
    return ax, fig

def plot_expt(labels=SELECTED_LABELS):
    expt_data = get_expt_energies()
    N_array = expt_data["N"].values
    energies = [expt_data[label].values for label in labels]
    ax, fig = plot_energies(N_array, energies)
    ax.set_title("Expt.")
    fig.savefig(SAVE_DIR / "expt_spectra.png")
    plt.close(fig)
    return None

def plot_calc(calc_data: dict[str, pd.DataFrame], min_n=MIN_N, max_n=MAX_N, labels=SELECTED_LABELS):
    N_array = np.arange(min_n, max_n+1, 2)
    for model_name, df in calc_data.items():
        energies = [df[label].values for label in labels]
        ax, fig = plot_energies(N_array, energies)
        ax.set_title(f"Model: {model_name}")
        (SAVE_DIR / model_name).mkdir(parents=True, exist_ok=True)
        try:
            fig.savefig(SAVE_DIR / f"{model_name}" / "spectra.png")
        finally:
            plt.close(fig)
    return None



if __name__ == "__main__":
    expt_data = get_expt_energies()
    calc_data: dict[str, pd.DataFrame] = get_calc_data()
    plot_expt()
    plot_calc(calc_data)