H_cc = 3.25         # -
C_S0 = 1            # mol/L (Initial concentration of sorbent)
z = 3               # m
R = 250e-6          # m (Radius of the fiber)   
#K_ext = 1.25e-3     # m/s (External mass transfer coefficient)

zeta = 250
D_CO2_l = 2e-10
u_l = D_CO2_l * z/(zeta*R**2)
print(f"Liquid velocity: {u_l:.5e} m/s")

Bi = 500
K_ext = Bi * H_cc * D_CO2_l / (R)  # m/s (External mass transfer coefficient)
print(f"External mass transfer coefficient: {K_ext:.5e} m/s")

theta = 0.24
D_S = theta * D_CO2_l  # m^2/s (Diffusion coefficient of salt in liquid phase)
print(f"Diffusion coefficient of sorbent in liquid phase: {D_S:.5e} m^2/s")
print(f'Diffusion coefficient of CO2 in liquid phase: {D_CO2_l:.5e} m^2/s')

C_CO2_g = 400*10**(-6)*101.325/(8.31446261815324*298)
epsilon = H_cc * C_CO2_g / C_S0
C_CO2_l = H_cc * C_CO2_g  # mol/l (Saturation concentration of CO2 in liquid phase)
print(f'CO2 concentration in gas phase: {C_CO2_g:.5e} mol/l')
print(f"Saturation concentration of CO2 in liquid phase: {C_CO2_l:.5e} mol/l")

Da = 30000
k_rxn = Da * D_CO2_l / (C_S0 * R**2)
print(f"Reaction rate constant: {k_rxn:.5e} m^3/(mol*s)")

