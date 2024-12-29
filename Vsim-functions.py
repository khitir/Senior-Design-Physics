                     #  D e f i n e       F u n c t i o n s
import numpy as np

def wz(z, w0, omega0):
    return np.sqrt(w0**2 * (1 + (z**2 / (w0**2 * omega0 / (2 * c))**2)))

def betaL(z, omega0, delta_omega, gamma, w0):
    return (z * gamma * delta_omega) / wz(z, w0, omega0)

def betaBAL(z, omega0, delta_omega, gamma, w0):
    return np.sqrt(1 + betaL(z, omega0, delta_omega, gamma, w0)**2)

def delta_omegaL(z, omega0, delta_omega, gamma, w0):
    return delta_omega / betaBAL(z, omega0, delta_omega, gamma, w0)

def omega0L(x, z, omega0, delta_omega, gamma, w0):
    return (omega0 + betaL(z, omega0, delta_omega, gamma, w0) / betaBAL(z, omega0, delta_omega, gamma, w0)**2 * x / (w0 * delta_omega))

def dOmega0L(x, z, omega0, delta_omega, gamma, w0):
    return omega0L(x, z, omega0, delta_omega, gamma, w0) - omega0

def Rz(z, w0, omega0):
    return z * (1 + ((omega0 * w0**2) / (2 * c))**2 / z**2)

def Rinvz(z, w0, omega0):
    return z / (z**2 + ((omega0 * w0**2) / (2 * c))**2)

def phi0L(x, y, z, omega0, delta_omega, gamma, w0, phi2In):
    omega0Loc = omega0L(x, z, omega0, delta_omega, gamma, w0)
    Rinv = Rinvz(z, w0, omega0)
    term1 = (Rinv * (y**2 + (x + z * gamma * (omega0 - omega0Loc))**2) * omega0Loc) / (2 * c)
    term2 = 0.5 * phi2In * (-omega0 + omega0Loc)**2
    term3 = (1 / c) * omega0Loc * (x * gamma * (-(omega0 - omega0Loc)) + z * (1 - 0.5 * gamma**2 * (-(omega0 - omega0Loc))**2))
    return term1 + term2 + term3



def phi1L(x, y, z, omega0, delta_omega, gamma, w0, phi2In):
    omega0Loc = omega0L(x, z, omega0, delta_omega, gamma, w0)
    Rinv = Rinvz(z, w0, omega0)
    term1 = -phi2In * (omega0 - omega0Loc)
    term2 = (1 / (2 * c)) * Rinv * (y**2 + (x + z * gamma * (omega0 - omega0Loc))**2 - 2 * z * gamma * (x + z * gamma * (omega0 - omega0Loc)) * omega0Loc)
    term3 = (1 / c) * ((x * gamma + z * gamma**2 * (omega0 - omega0Loc)) * omega0Loc + x * gamma * (-(omega0 - omega0Loc)) + z * (1 - 0.5 * gamma**2 * (-(omega0 - omega0Loc))**2))
    return term1 + term2 + term3

def phi2L(x, y, z, omega0, delta_omega, gamma, w0, phi2In):
    omega0Loc = omega0L(x, z, omega0, delta_omega, gamma, w0)
    Rinv = Rinvz(z, w0, omega0)
    term = (1 / (2 * c)) * (2 * x * gamma - 2 * Rinv * x * z * gamma + c * phi2In + 2 * z * gamma**2 * omega0 - 2 * Rinv * z**2 * gamma**2 * omega0 - 3 * z * gamma**2 * omega0Loc + 3 * Rinv * z**2 * gamma**2 * omega0Loc)
    return term




def At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In):
    omega0Loc = omega0L(x, z, omega0, delta_omega, gamma, w0)
    delta_omegaLoc = delta_omegaL(z, omega0, delta_omega, gamma, w0)
    phi0Loc = phi0L(x, y, z, omega0, delta_omega, gamma, w0, phi2In)
    phi1Loc = phi1L(x, y, z, omega0, delta_omega, gamma, w0, phi2In)
    phi2Loc = phi2L(x, y, z, omega0, delta_omega, gamma, w0, phi2In)
    betaBALoc = betaBAL(z, omega0, delta_omega, gamma, w0)
    w = wz(z, w0, omega0)
    zR = w0**2 * omega0 / (2 * c)
    phase_term = np.exp(-1j * omega0Loc * t)
    intensity_term = (w0 / w) * np.exp(-x**2 / (w**2 * betaBALoc**2)) * np.exp(-y**2 / w**2)
    phase_shift_term = np.exp(-1j * np.arctan(z / zR))
    phi_term = np.exp(1j * (phi0Loc - (delta_omegaLoc**2 * (t - phi1Loc)**2) / (4j + 4 * delta_omegaLoc**2 * phi2Loc)))
    normalization_factor = 2 * np.sqrt(np.pi) * np.sqrt(1 / delta_omegaLoc**2 - 1j * phi2Loc)
    return phase_term * intensity_term * phase_shift_term * phi_term / normalization_factor


def RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, Phi2In):

    w = wz(z, w0, omega0)
    zR = (omega0 * w0**2) / (2 * c)
    BetaBALoc = betaBAL(z, omega0, delta_omega, gamma, w0)
    Phi0Loc = phi0L(x, y, z, omega0, delta_omega, gamma, w0, Phi2In)
    Phi1Loc = phi1L(x, y, z, omega0, delta_omega, gamma, w0, Phi2In)
    Phi2Loc = phi2L(x, y, z, omega0, delta_omega, gamma, w0, Phi2In)
    omega0Loc = omega0L(x, z, omega0, delta_omega, gamma, w0)
    DeltaOmegaLoc = delta_omegaL(z, omega0, delta_omega, gamma, w0)
    
    term1 = np.cos(omega0Loc * t) * w0 / w
    term2 = np.exp(-(x**2 / (w**2 * BetaBALoc**2)))
    term3 = np.exp(-(y**2 / w**2))
    term4 = 1 / np.sqrt(1 + z**2 / zR**2)
    term5 = np.exp(-((DeltaOmegaLoc**2 * (t - Phi1Loc)**2) / (4 * (1 + DeltaOmegaLoc**4 * Phi2Loc**2))))
    term6 = np.cos(Phi0Loc - (DeltaOmegaLoc**4 * (t - Phi1Loc)**2 * Phi2Loc) / (4 * (1 + DeltaOmegaLoc**4 * Phi2Loc**2)))
    term7 = 1 / (2 * np.sqrt(np.pi) * np.sqrt(1 / DeltaOmegaLoc**2))
    
    result = term1 * term2 * term3 * term4 * (term5 * term6) * term7
    
    return result


###############################################################################################################################################

                     #       D e f i n e     u n i t s
# Time units
s = 1
ms = 10**-3 * s
us = 10**-6 * s
ns = 10**-9 * s
ps = 10**-12 * s
fs = 10**-15 * s

# Frequency units
Hz = 1
GHz = 10**9 * Hz
THz = 10**12 * Hz
MHz = 10**6 * Hz

# Length units
m = 1
cm = 10**-2 * m
mm = 10**-3 * m
um = 10**-6 * m
nm = 10**-9 * m
angs = 10**-10 * m
pm = 10**-12 * m

# Mass units
kg = 1
g = 10**-3 * kg

# Energy units
J = 1
mJ = 10**-3 * J
eV = 1 * J  # Assuming 'e' is the elementary charge
keV = 10**3 * eV
MeV = 10**6 * eV

# Power units
W = 1
kW = 10**3 * W
GW = 10**9 * W
mW = 10**-3 * W

# Pressure units
Torr = 1 / 760

# Temperature units
K = 1

# Constants
c = 299792458  # m/s
eps0 = 8.854187817 * 10**-12  # permittivity of free space
re = 2.81794092 * 10**-15  # m
e = 1.60217733 * 10**-19  # Coulomb, electron charge
Natm = 760 * 3.3 * 10**16 / cm**3  # no. particles in 1 atm
me = 9.1093897 * 10**-31  # kg
mp = me / 0.000544617013
hbar = 1.05457266 * 10**-34  # J s
h = 2 * 3.141592653589793 * hbar
kB = 1.3806503 * 10**-23  # J/K

#####################################################################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Check if each variable is defined and then delete it if so
if 'y' in locals():
    del y
if 'x' in locals():
    del x
if 'z' in locals():
    del z

# Define the parameters
lambda0 = 800*nm                   # wavelength
omega0 = 2 * np.pi * c / lambda0  #  frequency in rad/s
tau0 = 10*fs

delta_omega = 2 / tau0 # 1/fs 

theta_pf = np.pi/3

gamma = np.arctan(theta_pf) / omega0  # angular dispersion of focused beam

w0 = 1e-5  # Waist of the beam in meters
zc = 0  # Reference position in meters
phi2In = 0  # Input value for phi2

# Define the range of coordinates for the plot
x_range = np.linspace(-3*w0, 3*w0, 100)
y_range = np.linspace(-4*w0, 4*w0, 200)
z = 0  # Assuming z = 0 for simplicity
t = 0 

# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(y_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, y in enumerate(y_range):
        At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values))
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, y, z=0, t=0)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()


t = 100 * 1e-15   

# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(y_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, y in enumerate(y_range):
        At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values))
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, y, z=0, t=1e-13 )')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()


#######################################################################################################################################


# look at x-z plane

# Check if each variable is defined and then delete it if so
if 'y' in locals():
    del y
if 'x' in locals():
    del x
if 'z' in locals():
    del z

    
    
# Define the range of coordinates for the plot
x_range = np.linspace(-4*w0, 4*w0, 100)
z_range = np.linspace(-4*w0, 4*w0, 200)
t = 0  # Assuming t = 0 for simplicity

# Calculate the values of At for each (x, z) pair
At_values = np.zeros((len(z_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, z in enumerate(z_range):
        At_values[j, i] = At(x, 0, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values), extent=[x_range.min(), x_range.max(), z_range.min(), z_range.max()], aspect='auto')
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, y=0, z, t=0)')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.show()


t=100 * 1e-15
# Calculate the values of At for each (x, z) pair
At_values = np.zeros((len(z_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, z in enumerate(z_range):
        At_values[j, i] = At(x, 0, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values), extent=[x_range.min(), x_range.max(), z_range.min(), z_range.max()], aspect='auto')
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, y=0, z, t=1e-13)')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.show()


########################################################################################################################################################

# now we start using Real part of At()  function, compare if np.real(At)  and RealAt()  are same
import numpy as np
import matplotlib.pyplot as plt

# Check if each variable is defined and then delete it if so
if 'y' in locals():
    del y
if 'x' in locals():
    del x
if 'z' in locals():
    del z

# Define the parameters
lambda0 = 800*nm                   # wavelength
omega0 = 2 * np.pi * c / lambda0  #  frequency in rad/s
tau0 = 10*fs

delta_omega = 2 / tau0 # 1/fs 

theta_pf = np.pi/3

gamma = np.arctan(theta_pf) / omega0  # angular dispersion of focused beam

w0 = 1e-5  # Waist of the beam in meters
zc = 0  # Reference position in meters
phi2In = 0  # Input value for phi2

# Define the range of coordinates for the plot
x_range = np.linspace(-3*w0, 3*w0, 100)
y_range = np.linspace(-4*w0, 4*w0, 200)
z = 0  # Assuming z = 0 for simplicity
t = 0 #100 * 1e-15  

# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(y_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, y in enumerate(y_range):
        #At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        #At_values[j, i] = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        At_values[j, i] = np.real(At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In))

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values)) #extent=[x_range.min(), x_range.max(), y_range.min(), y_range.max()], aspect='auto')
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, y, z=0, t=0)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()


# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(y_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, y in enumerate(y_range):
        #At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        At_values[j, i] = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        #At_values[j, i] = np.real(At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In))

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values)) #extent=[x_range.min(), x_range.max(), y_range.min(), y_range.max()], aspect='auto')
plt.colorbar(label='|At|')
plt.title('My function Magnitude of At(x, y, z=0, t=0)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

##############################################################################################################################################


# Check if each variable is defined and then delete it if so
if 'y' in locals():
    del y
if 'x' in locals():
    del x
if 'z' in locals():
    del z


# Define the parameters
lambda0 = 800*nm                   # wavelength
omega0 = 2 * np.pi * c / lambda0  #  frequency in rad/s
tau0 = 10*fs

delta_omega = 2 / tau0 # 1/fs 

theta_pf = np.pi/3

gamma = np.arctan(theta_pf) / omega0  # angular dispersion of focused beam

w0 = 1e-5  # Waist of the beam in meters
zc = 0  # Reference position in meters
phi2In = 0  # Input value for phi2

# Define the range of coordinates for the plot
x_range = np.linspace(-3*w0, 3*w0, 100)
z_range = np.linspace(-4*w0, 4*w0, 200)
y = 0  # Assuming z = 0 for simplicity
t =0 #80 * 1e-15  

# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(z_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, z in enumerate(z_range):
        #At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        #At_values[j, i] = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        At_values[j, i] = np.real(At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In))

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values))
plt.colorbar(label='|At|')
plt.title('Magnitude of At(x, z, y=0, t=0)')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.show()


# Calculate the values of At for each (x, y) pair
At_values = np.zeros((len(z_range), len(x_range)), dtype=complex)
for i, x in enumerate(x_range):
    for j, z in enumerate(z_range):
        #At_values[j, i] = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        #At_values[j, i] = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
        At_values[j, i] = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Plot the result
plt.figure(figsize=(8, 6))
plt.imshow(np.abs(At_values))
plt.colorbar(label='|At|')
plt.title('My Function Magnitude of At(x, z, y=0, t=0)')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.show()


###########################################################################################################################################################
#  check values for same input, make sure they are same
# Check if each variable is defined and then delete it if so
if 'y' in locals():
    del y
if 'x' in locals():
    del x
if 'z' in locals():
    del z


####################################################################################################################################################
    
# Define input parameters
x = 0.1
y = -0.2
z = 0.3
t = 0.4
omega0 = 2 * np.pi * c / (800 * nm)  
tau0 = 10 * fs
delta_omega = 2 / tau0
theta_pf = np.pi / 3
gamma = np.arctan(theta_pf) / omega0
w0 = 1e-5
zc = 0
phi2In = 0

# Calculate outputs
At_result = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
RealAt_result = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)

# Compare results
real_part_At = np.real(At_result)
is_equal = np.allclose(real_part_At, RealAt_result)

# Print the comparison result
if is_equal:
    print("RealAt() is the real part of At().")
else:
    print("RealAt() is NOT the real part of At().")

####################################################################################################################################################
    


# Define ranges for x, y, and z
x_range = np.linspace(-3*w0, 3*w0, 10)  # Adjust the number of points as needed
y_range = np.linspace(-4*w0, 4*w0, 10)  # Adjust the number of points as needed
z_range = np.linspace(-3*w0, 3*w0, 10)         # Adjust the number of points as needed

# Initialize a flag to track the comparison result
all_equal = True

# Iterate over each combination of x, y, and z
for x in x_range:
    for y in y_range:
        for z in z_range:
            # Calculate outputs
            At_result = At(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
            RealAt_result = RealAt(x, y, z, t, omega0, delta_omega, gamma, w0, zc, phi2In)
            
            # Compare results
            real_part_At = np.real(At_result)
            is_equal = np.allclose(real_part_At, RealAt_result)
            #print(is_equal)
            # Update the flag if any comparison fails
            if not is_equal:
                all_equal = False
                break  # Exit the loop early if any comparison fails
                
        if not all_equal:
            break
    if not all_equal:
        break

# Print the overall comparison result
if all_equal:
    print("RealAt() is the real part of At() for all tested combinations of (x, y, z).")
else:
    print("RealAt() is NOT the real part of At() for some combinations of (x, y, z).")

