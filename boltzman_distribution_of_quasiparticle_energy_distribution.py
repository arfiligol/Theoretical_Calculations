import matplotlib.pyplot as plt
import numpy as np

# Constants
k_B = 1.38e-23  # Boltzmann constant in J/K
T_high = 300  # High temperature in K
T_low = 10  # Low temperature in K
E = np.linspace(0, 5e-20, 1000)  # Energy range in J

# Boltzmann distribution function
def boltzmann_distribution(E, T):
    return np.exp(-E / (k_B * T))

# Calculate distributions
P_high = boltzmann_distribution(E, T_high)
P_low = boltzmann_distribution(E, T_low)

# Normalize distributions
P_high /= np.trapz(P_high, E)
P_low /= np.trapz(P_low, E)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(E, P_high, label=f'High Temperature (T={T_high}K)', color='r')
plt.plot(E, P_low, label=f'Low Temperature (T={T_low}K)', color='b')
plt.xlabel('Energy (J)')
plt.ylabel('Probability Density')
plt.title('Boltzmann Energy Distribution at High and Low Temperatures')
plt.legend()
plt.grid(True)
plt.show()
