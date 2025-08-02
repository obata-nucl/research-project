import numpy as np
import pandas as pd
import torch
import os

from analyze import set_model
from config import BASE_DIR, MODEL_DIR
from NN_test import get_expt_data, run_with_timeout
from train import get_nuclei_range

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants
# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
# MIN_BETA : Minimum beta value for the nuclei
MIN_N = 86
MAX_N = 96
MIN_BETA = np.array([0.15, 0.20, 0.30, 0.35, 0.35, 0.35])
RMS_CRITERION = 2.0
HIDDEN_LAYERS = [128, 128, 128]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def Standardization(inputs_original: np.ndarray) -> np.ndarray:
    min = np.array([MIN_N, 2, -0.45])
    max = np.array([MAX_N, 7, 0.60])
    inputs_standardized = (inputs_original - min) / (max - min)
    return inputs_standardized

def prepare_test_data():
    Neutrons, Nu_bosons = get_nuclei_range(MIN_N, MAX_N)
    inputs_original = np.column_stack((Neutrons, Nu_bosons, MIN_BETA))
    inputs_scaled = Standardization(np.copy(inputs_original))
    return inputs_original, inputs_scaled

def save_results(save_path, N, Energies, ratio, params):
    # Check if the file exists, if not, create it and write the header
    if not os.path.exists(save_path):
        with open(save_path, "w") as f:
            f.write("N,2+_1,4+_1,6+_1,0+_2,E(4)/E(2),eps,kappa\n")
    # Append the data to the file
    with open(save_path, "a") as f:
        f.write(f"{int(N)},{Energies[0]},{Energies[1]},{Energies[2]},{Energies[3]},{ratio:.3f},{params[0]:.3f},{params[1]:.3f}\n")

def calc_RMS(th_energies: np.ndarray, expt_energies: np.ndarray) -> tuple:
    count = 0
    RMS = 0.0
    for i in range(len(th_energies)):
        if (not np.isnan(th_energies[i]) and not np.isnan(expt_energies[i])):
            count += 1
            RMS += (th_energies[i] - expt_energies[i]) ** 2
    return RMS, count

def evaluate(inputs_scaled, inputs_original, expt_data):
    model = set_model(HIDDEN_LAYERS)
    save_path = MODEL_DIR / "results_data" / "62" / f"{'_'.join(map(str, HIDDEN_LAYERS))}.csv"
    if os.path.exists(save_path):
        os.remove(save_path)

    total_count, energy_RMS = 0, 0.0
    ratio_count, ratio_RMS = 0, 0.0

    for i in range(len(inputs_scaled)):
        with torch.no_grad():
            input_tensor = torch.tensor(inputs_scaled[i], dtype=torch.float32).unsqueeze(0)
            output_tensor = model(input_tensor)
        sh_command = ["bash", str(MODEL_DIR / "src" / "NN_test.sh"), str(BASE_DIR / "NPBOS"), str(int(inputs_original[i, 0] + 62)), str(int(inputs_original[i, 1]))]
        sh_command.extend([f"{r:.3f}" for r in output_tensor.squeeze().numpy()])

        res = run_with_timeout(sh_command, timeout_sec=1.0)
        if res is None:
            print(f"Timeout expired for command: {' '.join(sh_command)}")
            continue
        try:
            energies = np.array([float(x) for x in res.split()])
        except ValueError:
            print(f"Error parsing output: {res}")
            continue

        RMS, count = calc_RMS(energies, expt_data[i])

        if count:
            energy_RMS += RMS
            total_count += count
        if energies[0] != 0 and not np.isnan(expt_data[i][0]):
            ratio_RMS += ((expt_data[i][1] / expt_data[i][0]) - (energies[1] / energies[0])) ** 2
            ratio_count += 1

            save_results(save_path, inputs_original[i, 0], energies, energies[1] / energies[0], output_tensor.squeeze().numpy())
        
    if total_count and ratio_count:
        total_RMS = np.sqrt(energy_RMS / total_count) + np.sqrt(ratio_RMS / ratio_count)
    if total_RMS < RMS_CRITERION:
        print(f"Model with hidden layers {HIDDEN_LAYERS} exceeds RMS criterion: {total_RMS:.3f}")

    return

if __name__ == "__main__":
    inputs_original, inputs_scaled = prepare_test_data()
    expt_N = get_expt_data()
    expt_data = expt_N.loc[expt_N["N"].isin(inputs_original[:, 0]), ["2+_1", "4+_1", "6+_1", "0+_2"]].to_numpy()
    evaluate(inputs_scaled, inputs_original, expt_data)