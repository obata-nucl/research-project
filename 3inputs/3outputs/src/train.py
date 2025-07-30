from model import NN, init_weights, loss_function

import torch
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import numpy as np
import itertools
import multiprocessing
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import pandas as pd
from torch.optim.lr_scheduler import ReduceLROnPlateau

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Constants

# INPUTS_DIM : Number of inputs
# OUTPUTS_DIM : Number of outputs
# NODES_PATTERN : Patterns of nodes in hidden layers
# HIDDEN_LAYERS : Number of hidden layers
# MODEL_DIR : Parent directory of the model
from config import BASE_DIR, MODEL_DIR, INPUTS_DIM, OUTPUTS_DIM, BATCH_SIZE, EPOCHS, LEARNING_RATE, WEIGHT_DECAY, NODES_PATTERN, HIDDEN_LAYERS, MIN_LR

MIN_N = 86
MAX_N = 96
VALIDATION_SIZE = 0.2
PATIENCE = 50 # For Early Stopping
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Settings
torch.set_num_threads(1)
# Set random seeds for reproducibility
def set_random_seeds(seed=42):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_random_seeds(42)


def generate_layer_configurations(nodes_pattern, layers_dims):
    all_combinations = set()
    for pattern in nodes_pattern:
        for layers_dim in layers_dims:
            all_combinations.update([tuple(c) for c in itertools.product(pattern, repeat=layers_dim)])
    return [list(c) for c in all_combinations]

def filter_layer_configurations(combinations):
    results = []
    for layers in combinations:
        is_valid = True
        for i in range(len(layers) - 1):
            current, next_ = layers[i], layers[i + 1]
            if max(current, next_) / min(current, next_) > 2:
                is_valid = False
                break
        if is_valid:
            results.append(layers)
    return results



# Scaling funcctions
def Standardization(inputs):
    return (inputs - np.min(inputs, axis=0)) / (np.max(inputs, axis=0) - np.min(inputs, axis=0))

def Normalization(inputs, mean=None, std=None):
    mean = np.mean(inputs, axis=0) if mean is None else mean
    std = np.std(inputs, axis=0) if std is None else std
    return (inputs - mean) / std



def get_nuclei_range(min_N, max_N) -> tuple[np.ndarray, np.ndarray]:
    # proton_numbers = np.arange(min_Z, max_Z + 1, 2)
    neutron_numbers = np.arange(min_N, max_N + 1, 2)
    magic_numbers = [2, 8, 20, 28, 50, 82, 126]
    # proton_bosons = []
    neutron_bosons = []
    for n in neutron_numbers:
        closest_magic_number = min(magic_numbers, key=lambda x : abs(n - x))
        neutron_bosons.append(abs(n - closest_magic_number) // 2)
    # for z in proton_numbers:
    #     closest_magic_number = min(magic_numbers, key=lambda x : abs(z - x))
    #     proton_bosons.append(abs(z - closest_magic_number) // 2)
    return neutron_numbers, np.array(neutron_bosons)

def prepare_train_data():
    Neutrons, Nu_bosons = get_nuclei_range(MIN_N, MAX_N)

    data_dir = BASE_DIR / "data" / "62"
    inputs, inputs_scaled, labels = [], [], []

    for N, nu_boson in zip(Neutrons, Nu_bosons):
        file_path = data_dir / f"{N}.csv"
        if not file_path.exists():
            print("train_data is not found :", file_path)
            break
        dataset = np.loadtxt(file_path, delimiter=",")
        beta = dataset[:, 0]
        energies = dataset[:, 1]
        energies -= np.min(energies)

        for i in range(len(beta)):
            inputs.append(np.array([N, nu_boson, beta[i]]))
            inputs_scaled.append(np.array([N, nu_boson, beta[i]]))
            labels.append(energies[i])
    
    X = torch.tensor(np.array(inputs), dtype=torch.float32)
    X_scaled = torch.tensor(Standardization(np.array(inputs_scaled)), dtype=torch.float32)
    Y = torch.tensor(np.array(labels), dtype=torch.float32)

    # Stratify by Neutron number (first column of X) to ensure balanced split
    stratify_labels = X[:, 0]
    X_train, X_val, X_scaled_train, X_scaled_val, Y_train, Y_val = train_test_split(
        X, X_scaled, Y, test_size=VALIDATION_SIZE, random_state=42, stratify=stratify_labels
    )

    return X_train, X_val, X_scaled_train, X_scaled_val, Y_train, Y_val



def train_worker(args):
    hidden_dims, X_train, X_val, X_scaled_train, X_scaled_val, Y_train, Y_val, process_id = args

    train_dataset = TensorDataset(X_scaled_train, X_train, Y_train)
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)
    
    val_dataset = TensorDataset(X_scaled_val, X_val, Y_val)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False, num_workers=0)

    model = NN(INPUTS_DIM, hidden_dims, OUTPUTS_DIM)
    model.apply(init_weights)  # Initialize weights
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    scheduler = ReduceLROnPlateau(optimizer, 'min', factor=0.1, patience=20)

    best_val_loss = float('inf')
    epochs_no_improve = 0
    best_model_state = None

    # --- Training loop ---
    for epoch in range(EPOCHS):
        # --- Training Phase ---
        model.train()
        for inputs_scaled, inputs_original, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs_scaled)
            loss = loss_function(outputs, inputs_original[:, 1], inputs_original[:, 2], labels)
            loss.backward()
            optimizer.step()
        
        # --- Validation Phase ---
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for inputs_scaled, inputs_original, labels in val_loader:
                outputs = model(inputs_scaled)
                loss = loss_function(outputs, inputs_original[:, 1], inputs_original[:, 2], labels)
                val_loss += loss.item()

        avg_val_loss = val_loss / len(val_loader)
        
        scheduler.step(avg_val_loss)

        # Check for improvement
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            epochs_no_improve = 0
            best_model_state = model.state_dict()
        else:
            epochs_no_improve += 1

        # if epochs_no_improve >= PATIENCE:
        #     break

        # current_lr = optimizer.param_groups[0]['lr']
        # if current_lr < MIN_LR:
        #     break
    
    print(f"[Process-{process_id}] Training finished. Config : {hidden_dims}, Best Val Loss : {best_val_loss:.4f}")

    if best_model_state:
        save_dir = MODEL_DIR / "NN_models"
        save_dir.mkdir(parents=True, exist_ok=True)
        save_path = save_dir / f"{'_'.join(map(str, hidden_dims))}.pth"
        torch.save(best_model_state, save_path)

    return {'config': hidden_dims, 'best_val_loss': best_val_loss}

def run_training():
    X_train, X_val, X_scaled_train, X_scaled_val, Y_train, Y_val = prepare_train_data()

    configs = generate_layer_configurations(NODES_PATTERN, HIDDEN_LAYERS)
    configs = filter_layer_configurations(configs)
    num_configs = len(configs)
    print(f"Total configurations: {num_configs}")

    tasks = [(config, X_train, X_val, X_scaled_train, X_scaled_val, Y_train, Y_val, i) for i, config in enumerate(configs)]

    total_cores = multiprocessing.cpu_count()
    num_processes = max(1, total_cores // 2)
    print(f"Starting training with {num_processes} processes...")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(train_worker, tasks)
    
    print("\nAll training processes completed.")

    # --- Save and display results ---
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='best_val_loss').reset_index(drop=True)
    
    results_dir = MODEL_DIR / "results_data"
    results_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(results_dir / "training_results.csv", index=False)

    print("\n--- Top 5 Configurations ---")
    print(results_df.head())
    print("\nResults saved to 'training_results.csv'")



if __name__ == "__main__":
    try:
        multiprocessing.set_start_method("fork", force=True)
    except RuntimeError:
        pass
    run_training()