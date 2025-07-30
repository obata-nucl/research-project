from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
MODEL_DIR = BASE_DIR / "2outputs"

# NN model parameters
INPUTS_DIM = 3
OUTPUTS_DIM = 2

# Training parameters
BATCH_SIZE = 64
EPOCHS = 500
# Optimizer parameters
LEARNING_RATE = 0.001
WEIGHT_DECAY = 1e-4
MIN_LR = 1e-6 # Minimum learning rate for early stopping

NODES_PATTERN = [[16, 32, 64], [32, 64, 128]]
HIDDEN_LAYERS = [3, 4, 5]
