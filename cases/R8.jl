# Case R8: vacuum | cold magnetized plasma in region VIII of the CMA diagram

freq = 5.8e9                                  # wave frequency [Hz]

medium_L = ColdPlasma()                       # z < 0: vacuum
medium_R = ColdPlasma(X = 8.0, Y = 4.0)       # z > 0: X = (ω_pe/ω)², Y = ω_ce/ω

nx = sind(45)                                 # tangential index nx = sin(θ_inc), θ_inc = 45°
ny = 0.0

A_in = [1.0, 1.0, 0.0, 0.0]                   # incoming amplitudes
