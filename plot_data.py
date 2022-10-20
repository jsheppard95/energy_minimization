#import python modules
import matplotlib.pyplot as plt
import pandas as pd

#NOTE:
#Everything below assumes unit atomic masses,
#such that forces = accelerations.

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
            # print(line)
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
                # print(key)
                if key not in result:
                    result[key] = [value]
                else:
                    result[key].append(value)
    return result


if __name__ == "__main__":
    # Read K = 100 energy minimization results file
    K100_file = "K100_energy_min.txt"
    leary_file = "LearyData.csv"
    result = ReadKEnergyMin(K100_file)

    # Plot minimum and average energies vs. N, for each value of K
    f, ax = plt.subplots()
    min_PE = ax.plot(result["N"], result["MinimumP.E"], ".", color="red",
                     label="Minimum")
    ax2 = ax.twinx()
    avg_PE = ax2.plot(result["N"], result["AverageP.E"], "x", color="blue",
                      label="Average")

    # Read and plot the expected data from Leary
    min_PE_exp = pd.read_csv(leary_file)
    #exp_PE = ax.plot(leary_file[])

    plots = min_PE + avg_PE
    labs = [p.get_label() for p in plots]
    ax.grid()
    ax.legend(plots, labs, loc=0)
    ax.set_title("K = 100")
    ax.set_xlabel(r"Number of Particles, $N$")
    ax.set_ylabel(r"Minimum Potential Energy (dimensionless), $U^*_{min}$")
    ax2.set_ylabel(r"Average Potential Energy (dimensionless), $U^*_{avg}$")
    f.savefig("MinAvgPE_N_K100.png")
    plt.show()
