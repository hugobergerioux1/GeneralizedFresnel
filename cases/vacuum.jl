# Case vacuum: vacuum on both sides (full transmission), to validate the waveports

freq = 5.8e9                                  # wave frequency [Hz]

medium_L = ColdPlasma()                       # z < 0: vacuum
medium_R = ColdPlasma()                       # z > 0: vacuum

nx = sind(45)                                 # tangential index nx = sin(θ_inc), θ_inc = 45°
ny = 0.0

A_in = [1.0, 1.0, 0.0, 0.0]                   # incoming amplitudes [A1 (TE), A2 (TM), A5, A6]
