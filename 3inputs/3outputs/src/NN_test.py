import numpy as np
import pandas as pd
import torch
import os
import subprocess
from model import NN
import signal

from config import BASE_DIR, MODEL_DIR, INPUTS_DIM, OUTPUTS_DIM, NODES_PATTERN, HIDDEN_LAYERS
from train import generate_layer_configurations, filter_layer_configurations, get_nuclei_range, Normalization

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants
# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
# MIN_BETA : Minimum beta value for the nuclei
MIN_N = 86
MAX_N = 96
MIN_BETA = np.array([0.15, 0.20, 0.30, 0.35, 0.35, 0.35])
RMS_CRITERION = 3.0
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def prepare_evaluate_data():
    Neutrons, Nu_bosons = get_nuclei_range(MIN_N, MAX_N)
    
    inputs_original = np.column_stack((Neutrons, Nu_bosons, MIN_BETA))
    inputs_scaled = Normalization(np.copy(inputs_original))
    return inputs_original, inputs_scaled

def get_expt_data():
    expt_path = BASE_DIR / 'data' / "62" / "expt.csv"
    df = pd.read_csv(expt_path)
    expt = df[(df["N"] >= MIN_N) & (df["N"] <= MAX_N)]
    return expt



def run_with_timeout(command, timeout_sec=1.0):
    try:
        proc = subprocess.Popen(
            list(map(str, command)),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            preexec_fn=os.setsid
        )
        stdout, _ = proc.communicate(timeout=timeout_sec)
        return stdout
    
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        return None
    
def calculate_rms(expt_data_N, energies):
    count = 0
    energy_RMS = 0
    for j in range(len(energies)):
        if np.isnan(expt_data_N[j]) or np.isnan(energies[j]):
            continue
        energy_RMS += (energies[j] - expt_data_N[j]) ** 2
        count += 1
    return energy_RMS, count
    


def evaluate_model(inputs_original, inputs_scaled, expt, configs):
    num_inputs = len(inputs_original)
    num_calc = len(configs)
    good_hidden_layers = []
    good_hidden_RMS = []


    for idx, config in enumerate(configs):
        hidden_RMS = 0
        ratio_RMS = 0
        total_count = 0
        ratio_count = 0
        model_path = MODEL_DIR / "NN_models" / f"{'_'.join(map(str, config))}.pth"
        if not model_path.exists():
            print(f"Model not found : {model_path}")
            continue

        model = NN(INPUTS_DIM, config, OUTPUTS_DIM)
        model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
        model.eval()

        for i in range(num_inputs):
            expt_data_N = expt.loc[expt["N"] == inputs_original[i, 0], ["2+_1", "4+_1", "6+_1", "0+_2"]].to_numpy().flatten()
            inputs_data = torch.tensor(inputs_scaled[i], dtype=torch.float32).unsqueeze(0)
            with torch.no_grad():
                outputs = model(inputs_data)
            
            results = [outputs[0][0], outputs[0][1], outputs[0][2]]
            sh_command = ["bash", str(MODEL_DIR / "src" / "NN_test.sh"), str(BASE_DIR / "NPBOS"), str(int(inputs_original[i, 0] + 62)), str(int(inputs_original[i, 1]))]
            sh_command.extend([f"{r:.3f}" for r in results])

            res = run_with_timeout(sh_command)
            if res is None:
                print(f"Timeout for command: {' '.join(sh_command)}")
                continue
            try:
                energies = list(map(float, res.split()))
            except ValueError as e:
                print(f"Error parsing output: {res} - {e}")
                continue

            energy_RMS, count = calculate_rms(expt_data_N, energies)
            if count:
                hidden_RMS += energy_RMS
                total_count += count
                if energies[0] != 0:
                    ratio_RMS += ((expt_data_N[1] / expt_data_N[0]) - (energies[1] / energies[0])) ** 2
                    ratio_count += 1
        hidden_RMS = np.sqrt(hidden_RMS / total_count) if total_count > 0 else 0
        ratio_RMS = np.sqrt(ratio_RMS / ratio_count) if ratio_count > 0 else 0

        if hidden_RMS + ratio_RMS < RMS_CRITERION:
            good_hidden_layers.append(config)
            good_hidden_RMS.append(hidden_RMS + ratio_RMS)
        if (idx + 1) % 10 == 0:
            print(f"Evaluated configuration {idx + 1}/{num_calc}")
    return good_hidden_layers, good_hidden_RMS

if __name__ == "__main__":
    inputs_original, inputs_scaled = prepare_evaluate_data()
    expt_energy_spectrum = get_expt_data()
    all_configs = generate_layer_configurations(NODES_PATTERN, HIDDEN_LAYERS)
    configs = filter_layer_configurations(all_configs)
    good_hidden_layers, good_hidden_RMS = evaluate_model(inputs_original, inputs_scaled, expt_energy_spectrum, configs)

    if not good_hidden_layers:
        print("No configurations met the RMS criterion. No valid models found.")
    else:
        sorted_results = sorted(zip(good_hidden_layers, good_hidden_RMS), key=lambda x: x[1])
        for layer, rms in sorted_results:
            print(f"Configuration: {layer}, RMS: {rms:.6f}")

        best_index = np.argmin(good_hidden_RMS)
        print(f"Best configuration: {good_hidden_layers[best_index]}, RMS: {good_hidden_RMS[best_index]:.6f}")