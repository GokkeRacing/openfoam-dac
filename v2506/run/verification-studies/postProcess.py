import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def get_second_highest_numbered_folder(path='.'):
    p = Path(path)
    numeric_folders = []
    for f in p.iterdir():
        if f.is_dir():
            try:
                val = int(f.name)
            except ValueError:
                try:
                    val = float(f.name)
                except ValueError:
                    continue
            numeric_folders.append((val, f.name))
    if len(numeric_folders) < 2:
        return None
    numeric_folders.sort(key=lambda x: x[0], reverse=True)
    return numeric_folders[1][1]

def read_openfoam_scalar_field(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the start of the nonuniform List<scalar>
    for i, line in enumerate(lines):
        if 'nonuniform List<scalar>' in line:
            n_values = int(lines[i+1].strip())
            start_idx = i + 3  # Skip the '(' line
            break
    else:
        raise ValueError("Could not find nonuniform List<scalar> in file.")

    # Read the values
    values = []
    for line in lines[start_idx:start_idx + n_values]:
        values.append(float(line.strip()))

    return np.array(values)


if __name__ == "__main__":
    cwd = Path(__file__).parent
    folder_latestTime = cwd /get_second_highest_numbered_folder(cwd)

    data_CO2 = read_openfoam_scalar_field(folder_latestTime / 'CO2')
    data_S = read_openfoam_scalar_field(folder_latestTime / 'S')
    data_V = read_openfoam_scalar_field(folder_latestTime / 'V')
    data_Cy = read_openfoam_scalar_field(folder_latestTime / 'Cy')

    data = np.column_stack((data_Cy, data_CO2, data_S, data_V))  # add more fields as needed

    # Get unique Cy values and the indices to reconstruct the groups
    unique_Cy, inverse_indices = np.unique(data_Cy, return_inverse=True)

    # Prepare an array to hold the means
    means = np.zeros((len(unique_Cy), data.shape[1]))

    for i, cy_val in enumerate(unique_Cy):
        group = data[inverse_indices == i]
        weights = group[:, 3]  # V is the 4th column
        means[i, 0] = np.average(group[:, 0], weights=weights)  # Cy
        means[i, 1] = np.average(group[:, 1], weights=weights)  # CO2
        means[i, 2] = np.average(group[:, 2], weights=weights)  # S
        means[i, 3] = np.average(group[:, 3], weights=weights)  # V

    paper_phi_l = np.loadtxt('phi_l.txt')
    paper_phi_s = np.loadtxt('phi_s.txt')

    plt.figure(0)
    plt.plot(means[:, 0]/3*250, means[:, 1]/np.max(means[:, 1]), linewidth=2, color='tab:red', label='Simulation')
    plt.plot(paper_phi_l[:, 0], paper_phi_l[:, 1], 's', linewidth=2, color='tab:red',label='Diederichsen and Hatton (2022)')
    plt.legend(loc="upper left",fontsize=12)
    plt.xlabel(r'Non-dimensional axial position, $\zeta$',fontsize=12)
    plt.ylabel(r'Non-dimensional radial-averaged concentration, $\phi_{\text{CO}_2}$',fontsize=12)
    plt.grid()
    plt.savefig("CO2.svg")

    plt.figure(1)
    plt.plot(means[:, 0]/3*250, means[:, 2]/np.max(means[:, 2]), linewidth=2, color='tab:blue', label='Simulation')
    plt.plot(paper_phi_s[:, 0], paper_phi_s[:, 1], 's', linewidth=2, color='tab:blue',label='Diederichsen and Hatton (2022)')
    plt.legend(loc="upper right",fontsize=12)
    plt.xlabel(r'Non-dimensional axial position, $\zeta$',fontsize=12)
    plt.ylabel(r'Non-dimensional radial-averaged concentration, $\phi_\text{S}$',fontsize=12)
    plt.grid()
    plt.savefig("S.svg")
