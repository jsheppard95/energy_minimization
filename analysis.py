#import python modules
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit

def ReadKEnergyMin(file):
    """
    Reads Conjugate Gradient output file and returns N values and minimum,
    average, and maximum potential energies.

    Parameters:
    -----------
    file : str ; path to energy minimization over K file

    Returns:
    --------
    result : dict ; dictionary of N, minimum P.E, and maximum P.E for the data
        in `file`
        {"N": [N values], "MinimumP.E": [MinPEVals],
         "AverageP.E": [AvgPEVals], "MaximumP.E": [MaxPEVals]}
    """
    # Data to save
    result = {}
    # Open file for reading
    with open(file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            entries = line.split(",")
            # Result
            # ['N: 26', ' Minimum P.E: -101., ' Average P.E: -6.',
            #  ' Maximum P.E: -0.7\n']
            for entry in entries:
                data = entry.split(":")
                key = data[0].replace(" ", "")
                if key == "N":
                    value = int(data[1])
                else:
                    value = float(data[1])
                if key not in result:
                    result[key] = [value]
                else:
                    result[key].append(value)
    return result


def PlotResults(K_energy_min_file, K):
    # Read K energy minimization results file
    result = ReadKEnergyMin(K_energy_min_file)


    # Plot minimum and average energies vs. N, for each value of K
    f, ax = plt.subplots()
    min_PE = ax.plot(result["N"], result["MinimumP.E"], ".", color="black",
                     label="Minimum")
    ax2 = ax.twinx()
    avg_PE = ax2.plot(result["N"], result["AverageP.E"], "x", color="blue",
                      label="Average")
    
    # Plot expected data
    leary_fname = "LearyData.csv"
    leary_file = os.path.join("data", leary_fname)
    leary_df = pd.read_csv(leary_file)
    leary_PE = ax.plot(leary_df["N"], -leary_df["minusUmin"], "^", color="black",
                     label="Leary Data")

    # Multi-axis plot customization
    plots = min_PE + leary_PE + avg_PE
    labs = [p.get_label() for p in plots]
    ax.grid()
    ax.legend(plots, labs, loc=0)
    ax.set_title(f"K = {K}")
    ax.set_xlabel(r"Number of Particles, $N$")
    ax.set_ylabel(r"Minimum/Leary Potential Energy (dimensionless), $U^*_{min}$")
    ax2.set_ylabel(r"Average Potential Energy (dimensionless), $U^*_{avg}$")
    ax2.yaxis.label.set_color("blue")
    output_dir = "plots"
    outfile_name = os.path.basename(K_energy_min_file).replace(".txt", ".png")
    f.savefig(os.path.join(output_dir, outfile_name))


def UMacro(N, a, b, c):
    """
    Model for U_macro:
    U_macro = a + bN^(2/3) + cN
    Cluster Energy ~ Surface Area (Surface Tension) and Volume (bulk energy
    density)
    
    Parameters:
    -----------
    N : int ; number of particles
    a : float ; fit parameter
    b : float ; fit parameter
    c : float ; fit parameter

    Returns:
    --------
    U_macro : float ; global energy minimum 
    """
    return a + b*N**(2/3) + c*N

def FitUMacro(K_energy_min_file, K):
    pass


if __name__ == "__main__":
    # Set file paths
    K100_fname = "K100_energy_min.txt"
    K100_file = os.path.join("data", K100_fname)
    # Plot K100 results
    PlotResults(K100_file, 100)

    # Fit Umacro
    K100_data = ReadKEnergyMin(K100_file)
    popt, pcov = curve_fit(UMacro, K100_data["N"], K100_data["MinimumP.E"])

    # Plot result (fit and actual N)
    N_vals = np.linspace(2, 26, 1000)
    U_pred = UMacro(N_vals, *popt)
    f1, ax1 = plt.subplots()
    act = ax1.plot(K100_data["N"], K100_data["MinimumP.E"], ".",
                  label=r"Actual, $U^*_{min}$")
    pred = ax1.plot(N_vals, U_pred, label=r"Predicted $U_{macro}$")

    # Plot customization
    def style_plot(ax, K, param_str):
        """
        Helper function to do common styling
        """
        ax.grid()
        ax.set_xlabel(r"Number of Particles, $N$")
        ax.set_title(f"K = {K}\n"
                     r"$U_{macro} = a + bN^{2/3} + cN$")
        anchored_text = AnchoredText(param_str, loc="lower left")
        ax.add_artist(anchored_text)

    plots = act + pred
    labs = [p.get_label() for p in plots]
    param_str = f"a = {popt[0]:.2f}, b={popt[1]:.2f}, c={popt[2]:.2f}"
    ax1.legend(plots, labs, loc=0)
    ax1.set_ylabel(r"Minimum Potential Energy (dimensionless), $U^*_{min}$")
    style_plot(ax1, 100, param_str)

    # Plot difference U_min - U_macro
    N_disc = np.arange(2, 27)
    U_pred_disc = UMacro(N_disc, *popt)
    resid = K100_data["MinimumP.E"] - U_pred_disc
    f2, ax2 = plt.subplots()
    ax2.plot(N_disc, resid, ".")
    ax2.set_ylabel(r"Residuals, $U^*_{min} - U_{macro}$")
    style_plot(ax2, 100, param_str)

    # Save plots
    output_dir = "plots"
    fit_outfile_name = os.path.basename(K100_file).replace(".txt", "_UmacroFit.png")
    f1.savefig(os.path.join(output_dir, fit_outfile_name))
    resid_outfile_name = os.path.basename(K100_file).replace(".txt", "_UmacroResid.png")
    f2.savefig(os.path.join(output_dir, resid_outfile_name))
    plt.show()
