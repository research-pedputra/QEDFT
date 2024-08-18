from ase.io.cube import read_cube_data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class ChargeDensityPlotter:
    def __init__(self, cube_files, dimensions, scale_factor):
        self.cube_files = cube_files
        self.dimensions = dimensions
        self.scale_factor = scale_factor
        self.data = []
        self.zq_CDNSDEC = None

    def load_cube_data(self):
        """Load data from the provided .cube files."""
        for file in self.cube_files:
            data, _ = read_cube_data(file)
            self.data.append(data)
    
    def compute_charge_density_difference(self):
        """Compute the charge density difference."""
        if len(self.data) < 3:
            raise ValueError("At least three data sets are required.")
        diff = (self.data[0] - self.data[1] - self.data[2])
        volume_element = np.prod(self.dimensions)
        self.CDNSDEC = (diff / volume_element) * self.scale_factor

    def reshape_and_sum_charge_density(self):
        """Reshape and sum the charge density along the z-axis."""
        q_CDNSDEC = np.reshape(self.CDNSDEC, self.dimensions)
        self.zq_CDNSDEC = [np.sum(q_CDNSDEC[:, :, int(i)]) for i in range(self.dimensions[2])]

    def plot_charge_density_difference(self, output_file):
        """Plot the charge density difference and save the plot."""
        if self.zq_CDNSDEC is None:
            raise ValueError("Charge density data has not been computed yet.")
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.plot(self.zq_CDNSDEC, color='red')
        ax.set_xlabel('z-distance ($\AA$)', fontsize=16, fontweight='bold')
        ax.set_ylabel('Total Charge Density Difference (e/$\AA$)', fontsize=16, fontweight='bold', color='red')
        ax.tick_params(axis='x', which='major', labelsize=16)
        ax.tick_params(axis='y', which='major', labelsize=16)
        ax.axhline(y=0, linewidth=0.5, color='k', linestyle=(0, (8, 10)))
        ax.set_xticklabels(np.linspace(0, self.dimensions[2] * 0.127, 5))  # Approximate angstrom scale
        ax.set_xticks(np.linspace(0, self.dimensions[2], 5))
        ax.spines['left'].set_color('red') 
        ax.tick_params(axis='y', colors='red')
        ax.set_ylim(-0.08, 0.08)
        plt.savefig(output_file, dpi=300, transparent=True)
        plt.show()

def main():
    plt.rcParams.update({'font.family':'Times New Roman'})
    plt.rcParams['figure.figsize'] = [8, 6]

    # Parameters passed separately
    cube_files = [
        'D:/Research/Charge_Density_Cal_110_SnO2_NiO_CO2_S5/All/110_SnO2_NiO+CO2.cube',
        'D:/Research/Charge_Density_Cal_110_SnO2_NiO_CO2_S5/All/110_SnO2_NiO.cube',
        'D:/Research/Charge_Density_Cal_110_SnO2_NiO_CO2_S5/All/110_CO2.cube'
    ]
    dimensions = (160, 150, 200)  # FFT grid dimensions
    scale_factor = 20.0973892212 * 19.1100006104 * 25.4062004089 * (0.5291772108 ** 3)

    # Create the plotter instance with parameters
    plotter = ChargeDensityPlotter(cube_files, dimensions, scale_factor)
    plotter.load_cube_data()
    plotter.compute_charge_density_difference()
    plotter.reshape_and_sum_charge_density()
    plotter.plot_charge_density_difference('Charge_transfer_110_SnO2+NiO_S5.png')

if __name__ == "__main__":
    main()
