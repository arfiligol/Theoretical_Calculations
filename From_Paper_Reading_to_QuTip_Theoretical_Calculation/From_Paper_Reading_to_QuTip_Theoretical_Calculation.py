from qutip import Qobj, basis, num, qzero, destroy,identity
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, hbar, giga

# Define the parameters value
EJ = 17
EC = 1 
P = -1  # For simplicity, consider P = -1 (odd parity)
omega_r = 9.202 * 2 * np.pi * giga  # Resonator frequency in rad/s
g = 0.01 * 2 * np.pi* giga  # Coupling strength in Hz, adjust as needed
num_levels_to_observe = 5  # Number of energy levels to observe

# Define the number of states (N+1 states, from 0 to N)
N = 4 # Totally 11 states (include zero state).

# Create basis states |n> for n = 0 to N
basis_states = [basis(N+1, n) for n in range(N+1)]


N_hat= destroy(N+1).dag()*destroy(N+1)



# Initialize the \hat{\varphi} operator as a zero matrix of size (N+1)x(N+1)
phi_hat = qzero(N + 1)

# Construct the \hat{\varphi} operator
for m in range(N+1):
    if m + 1 <= N: 
        phi_hat += (1j / 2) * basis_states[m] * basis_states[m + 1].dag()
    if m - 1 >= 0:
        phi_hat += (-1j / 2) * basis_states[m] * basis_states[m - 1].dag()

# Initialize the shift operators
e_iphi = qzero(N + 1)
e_imphi = qzero(N + 1)

# Construct the shift operator e^{i \hat{\varphi}}
for n in range(N):
    e_iphi += basis(N + 1, n + 1) * basis(N + 1, n).dag()

# Construct the shift operator e^{-i \hat{\varphi}}
for n in range(1, N + 1):
    e_imphi += basis(N + 1, n - 1) * basis(N + 1, n).dag()

# Construct cos(\hat{\varphi})
cos_phi = 0.5 * (e_iphi + e_imphi)

# Constructing the Hamiltonian H_{CPB}
# For different n_g, we will have different H_CPB, so convert it to a function which return different H_CPB

def return_H_CPB(n_g, P) -> Qobj:    
    return (4 * EC * (N_hat - n_g + (P - 1)/4)**2) - (EJ * cos_phi)

def hamiltonian(ng,P):

    m = np.diag(4 * EC * (np.arange(-N,N+1)-ng+ (P - 1)/4)**2) + 0.5 * EJ * (np.diag(-np.ones(2*N), 1) + 
                                                               np.diag(-np.ones(2*N), -1))
    return Qobj(m)

#%%
# Initialize lists to collect eigenenergies for different P
eigenenergies_list_p_neg1 = []
eigenenergies_list_p_pos1 = []

# Define the range of n_g values
n_g_values = np.linspace(0, 1, 101)

# Compute eigenenergies for each n_g and P = -1
for n_g in n_g_values:
    H_CPB = hamiltonian(n_g,-1)#return_H_CPB(n_g, -1)
    eigenenergies = H_CPB.eigenenergies()
    eigenenergies_list_p_neg1.append(eigenenergies / EJ)

# Compute eigenenergies for each n_g and P = 1
for n_g in n_g_values:
    H_CPB = hamiltonian(n_g,1)#return_H_CPB(n_g, 1)
    eigenenergies = H_CPB.eigenenergies()
    eigenenergies_list_p_pos1.append(eigenenergies / EJ)

# Convert the lists to numpy arrays for easier plotting
eigenenergies_array_p_neg1 = np.array(eigenenergies_list_p_neg1)
eigenenergies_array_p_pos1 = np.array(eigenenergies_list_p_pos1)

# Plot the eigenenergies as a function of n_g for P = -1 (solid lines)
for i in range(num_levels_to_observe):
    plt.plot(n_g_values, eigenenergies_array_p_neg1[:, i], label=f'Level {i} (P=-1)', linestyle='-')

# Plot the eigenenergies as a function of n_g for P = 1 (dashed lines)
for i in range(num_levels_to_observe):
    plt.plot(n_g_values, eigenenergies_array_p_pos1[:, i], label=f'Level {i} (P=1)', linestyle='--')

plt.xlabel('$n_g$')
plt.ylabel('Eigenenergies / $E_J$')
plt.legend()
plt.show()
