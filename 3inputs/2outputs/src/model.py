import torch
import torch.nn as nn
import torch.nn.functional as F

# NNのモデルを定義
class NN(nn.Module):
    def __init__(self, input_dim, hidden_dims, output_dim, dropout_rate=0.2):
        super(NN, self).__init__()
        layers = []
        preveous_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(preveous_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout_rate))
            preveous_dim = hidden_dim
        layers.append(nn.Linear(preveous_dim, output_dim))
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        output = self.model(x)
        eps = F.softplus(output[:, 0])
        kappa = -F.softplus(output[:, 1])
        return torch.stack([eps, kappa], dim=1)

def init_weights(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
            nn.init.zeros_(m.bias)
    # elif isinstance(m, nn.Conv2d):
    #     nn.init.kaiming_uniform_(m.weight, nonlinearity='relu')
    #     if m.bias is not None:
    #         nn.init.zeros_(m.bias)
    # elif isinstance(m, nn.BatchNorm2d):
    #     nn.init.ones_(m.weight)
    #     nn.init.zeros_(m.bias)



def IBM2_energy(output, n_nu, beta_f):
    # IBM parameters
    eps = output[:, 0]
    kappa = output[:, 1]

    beta_b = beta_f * 5
    beta2 = beta_b * beta_b
    num_bosons = 6 + n_nu
    deno = 1 / (1 + beta2)

    # Calculate the IBM2 energy
    IBM2_energy = (
        eps * num_bosons * beta2 * deno
    ) + (
        6 * n_nu * kappa * beta2 * deno * deno * (4 + 1.3 * 2 * torch.sqrt(torch.tensor(2 / 7)) * beta_b + 0.8 * beta2 / 7)
    )

    return IBM2_energy

def loss_function(output, n_nu, beta_f, HFB_energy):
    predicted_energy = IBM2_energy(output, n_nu, beta_f)
    loss = F.mse_loss(predicted_energy, HFB_energy)
    return loss
