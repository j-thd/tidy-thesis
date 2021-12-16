import os
import generate_optimization_results as gor
import aggregrate_optimization_results as aor
import optimization.settings

# Optimization and model settings
settings = optimization.settings.settings_1D_rectangular_multichannel # All settings, such as model resolution, Nusselt Relations, pressure drop relations etc..
F_desired = 1e-3 # [N] Desired thrust. Earlier results may be overwritten if run at same thrust level
# Settings for where to store results. Between pre- and suffix the thrust level is added.
results_folder_prefix = "optimization_results/" # Results for one thrust level will be stored in a folder beginning with this name
results_folder_suffix = "-alpha/" # A suffix, which may be convenient if multiple settings were used, and previous results must be kept
if not results_folder_suffix.endswith('/'):
    print("Results folder suffix must end with '/'")
    exit()


# Range of chamber temperatures to evaluate (multiples of 25K are nice to work with, to get nice round chamber temperatures)
T_chamber_min = 500  # [K] Lowest temperature to evaluate (must be higher than saturation temperature for pressure specified in settings)
T_chamber_max = 1100 # [K] Highest tempereature to evaluate
delta_T = 25 # [K] Step size between temperatures. Preferably, above range is a multiple of this value
steps = round((T_chamber_max-T_chamber_min)/delta_T)+1 # Amount of steps


results_folder = results_folder_prefix + "{:3.5f}mN".format(F_desired*1e3) + results_folder_suffix # [str] The final complete name of the results folder
# Check if the folder exists. Halt if it does.
skip_generating_results = False
if os.path.isdir(results_folder):
    print("The following folder already exists: {}".format(results_folder))
    print("Delete the folder manually if you are sure about overwriting results")
    print("Continuing to plotting previous results")
    skip_generating_results = True
else:
    os.mkdir(path=results_folder)
    print("Storing results in: {}".format(results_folder))


# Range of number of channels to evaluate
# NOTE: Main reason to have this available is to reduce the range for computation speed,
#   but the optimization can also get stuck at lower thrust levels and too high channel counts (due to numerical issues)
number_of_channels_low = 1 # [-] Lowest number of channels
number_of_channels_high = 10 # [-] Highest number of channels to evaluate

if not skip_generating_results:
    print("Generating results for F = {:1.5f} mN".format(F_desired*1e3))
    gor.run(
        F_desired=F_desired,
        T_low=T_chamber_min,
        T_high=T_chamber_max,
        channels_min = number_of_channels_low,
        channels_max = number_of_channels_high, 
        steps=steps,
        str_folder=results_folder,
        settings=settings)

## Plotting results
aor.plot_results(str_folder=results_folder)