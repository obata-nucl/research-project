import numpy as np
import torch
import os

from model import NN
from config import BASE_DIR, MODEL_DIR, INPUTS_DIM, OUTPUTS_DIM
from NN_models import prepare_evaluate_data, run_with_timeout

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants
# MIN_N : Minimum of the neutron number in nuclei
# MAX_N : Maximum of the neutron number in nuclei
# MIN_BETA : Minimum beta value for the nuclei
MIN_N = 86
MAX_N = 96
MIN_BETA = np.array([0.15, 0.20, 0.30, 0.35, 0.35, 0.35])

OUTPUT_LABELS = ["eps", "kappa", "chi_nu"]

GOOD_MODELS = [[16, 16, 16, 16, 16]
               ]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def set_model(hidden_dim):
    model_path = MODEL_DIR / "NN_models" / f"{'_'.join(map(str, hidden_dim))}.pth"
    if not model_path.exists():
        print(f"Model not found: {model_path}")
        return None
    model_state = torch.load(model_path, map_location="cpu")
    model = NN(INPUTS_DIM, hidden_dim, OUTPUTS_DIM)
    model.load_state_dict(model_state)
    return model

def save_results(file_path, N, Energies, ratio, params):
    # Check if the file exists, if not, create it and write the header
    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            f.write("N,2+_1,4+_1,6+_1,0+_2,E(4)/E(2),eps,kappa,chi_nu\n")
    # Append the data to the file
    with open(file_path, "a") as f:
        f.write(f"{int(N)},{Energies[0]},{Energies[1]},{Energies[2]},{Energies[3]},{ratio:.3f},{params[0]:.3f},{params[1]:.3f},{params[2]:.3f}\n")

def analyze_top5(good_models, inputs_original, inputs_scaled):
    for hidden_dim in good_models:
        print(f"Hidden dims: {hidden_dim}")
        model = set_model(hidden_dim)
        if model is None:
            continue
        save_path = MODEL_DIR / "results_data" / "62" / f"{'_'.join(map(str, hidden_dim))}.csv"

        # Delete the file if it exists
        if os.path.exists(save_path):
            os.remove(save_path)

        for i in range(len(inputs_scaled)):
            with torch.no_grad():
                output = model(torch.tensor(inputs_scaled[i], dtype=torch.float32).unsqueeze(0))
            results = [output[0][0].item(), output[0][1].item(), output[0][2].item()]
            command = [
                "bash", str(MODEL_DIR / "src" / "NN_test.sh"),
                str(BASE_DIR / "NPBOS"), str(int(inputs_original[i, 0] + 62)),
                str(int(inputs_original[i, 1]))
            ]
            for result in results:
                command.append(f"{result:.3f}")
            res = run_with_timeout(command, timeout_sec=1.0)
            if res is None:
                print(f"Timeout for command: {' '.join(command)}")
                print("Exiting program due to timeout.")
                exit(1)  # Exit the program on timeout
            try:
                energies = list(map(float, res.split()))
                ratio = energies[1] / energies[0] if energies[0] != 0 else float('inf')
                print(energies, ratio)
                save_results(save_path, inputs_original[i, 0], energies, ratio, results)
            except ValueError as e:
                print("Failed tp pause output:", e)
                continue
        print("======================================")
    return None

if __name__ == "__main__":
    inputs_original, inputs_scaled = prepare_evaluate_data()
    analyze_top5(GOOD_MODELS, inputs_original, inputs_scaled)