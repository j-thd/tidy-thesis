# File to automatically run optimization for various temperatures
import main_1D_rectangular_multichannel_optimization as optim
import sys
import numpy as np

def run(F_desired, T_low, T_high, channels_min, channels_max, settings, steps, str_folder):
    T_range = np.linspace(start=T_low, stop=T_high, num=steps)
    T_iter = np.nditer(T_range)
    for T in T_iter:
        optim.run(F_desired=F_desired, T_chamber=round(float(T)), channels_min=channels_min, channels_max=channels_max, settings=settings, str_folder=str_folder)

if __name__ == "__main__":
    run()