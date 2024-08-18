import numpy as np
import matplotlib.pyplot as plt

class BandEnergy:
    def __init__(self):
        self.band_data = None
        self.fermi_energy = None
        self.k = None
        self.bands = None
        self.high_symmetry_points = None
        self.fig, self.ax = None, None

    def get_bandstructure(self, data, fermi_energy, kpoints):
        self.band_data = data
        self.fermi_energy = fermi_energy
        self.k = np.unique(data[:, 0])
        self.bands = np.reshape(data[:, 1], (-1, len(self.k)))
        self.high_symmetry_points = kpoints

    def plot_band_structure(self):
        if self.band_data is None:
            raise ValueError("Band structure data has not been loaded.")
        self.fig, self.ax = plt.subplots(1, 1, figsize=(6, 6))
        for i, band in enumerate(self.bands):
            self.ax.plot(self.k, band - self.fermi_energy, linewidth=1, alpha=0.5, color='r')
        self.ax.set_xlim(min(self.k), max(self.k))
        self.ax.set_ylim(-5, 5)
        
        if self.high_symmetry_points is not None:
            x_coords = [x_coord for _, x_coord in self.high_symmetry_points]
            x_labels = [symmetry_point for symmetry_point, _ in self.high_symmetry_points]
            self.ax.set_xticks(x_coords)
            self.ax.set_xticklabels(x_labels)

            for _, x_coordinate in self.high_symmetry_points:
                self.ax.axvline(x=x_coordinate, color='k', linestyle='-', linewidth=1)
        self.ax.set_ylabel("E-E$_f$ (eV)", fontweight='bold')

    def plot_band_number(self, band_number, overlay=True):
        if band_number < 1 or band_number > len(self.bands):
            raise ValueError("Invalid band number.")
        
        if overlay:
            if self.fig is None or self.ax is None:
                self.plot_band_structure()  # Plot the band structure if not plotted already
            self.ax.plot(self.k, self.bands[band_number - 1] - self.fermi_energy, linewidth=1, color='b')
        else:
            fig, ax = plt.subplots(1, 1, figsize=(6, 6))
            for i, band in enumerate(self.bands):
                ax.plot(self.k, band - self.fermi_energy, linewidth=1, alpha=0.5, color='r')
            ax.plot(self.k, self.bands[band_number - 1] - self.fermi_energy, linewidth=1, color='b')
            ax.set_xlim(min(self.k), max(self.k))
            ax.set_ylim(-5, 5)
            
            if self.high_symmetry_points is not None:
                x_coords = [x_coord for _, x_coord in self.high_symmetry_points]
                x_labels = [symmetry_point for symmetry_point, _ in self.high_symmetry_points]
                ax.set_xticks(x_coords)
                ax.set_xticklabels(x_labels)
                for _, x_coordinate in self.high_symmetry_points:
                    ax.axvline(x=x_coordinate, color='k', linestyle='-', linewidth=1)
            ax.set_ylabel("E-E$_f$ (eV)", fontweight='bold')
            ax.set_title(f"Band {band_number}")
            plt.show()
    
    def shift_band(self, start_band_number, shift_value):
        if start_band_number < 1 or start_band_number > len(self.bands):
            raise ValueError("Invalid band number.")
        # Shift the specified band and all subsequent bands
        for i in range(start_band_number - 1, len(self.bands)):
            self.bands[i] += shift_value
    
    def plot_shifted_band_structure(self, start_band_number, shift_value):
        # Apply the shift
        self.shift_band(start_band_number, shift_value)
        # Plot the entire band structure with the applied shift
        self.plot_band_structure()
    
    def get_band_data_point(self, band_number):
        if band_number < 1 or band_number > len(self.bands):
            raise ValueError("Invalid band number.")
        
        return self.k, self.bands[band_number - 1] - self.fermi_energy
    
    def load_pdos_data(self, file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()

        data_start_index = None
        for i, line in enumerate(lines):
            if line.startswith("#"):
                data_start_index = i + 1
                break

        if data_start_index is None:
            raise ValueError("Invalid file format. Data not found after header.")

        data = np.loadtxt(file_path, skiprows=data_start_index)
        return data[:, 0], data[:, 1:]

    def get_aopdos(self, files, fermi_energy):
        aopdos_data = {}
        for name, (color, file_path) in files.items():
            energy, pdos = self.load_pdos_data(file_path)
            shifted_energy = energy - fermi_energy
            aopdos_data[name] = {"Energy": shifted_energy, "PDOS": pdos, "Color": color}
        return aopdos_data

    def plot_aopdos(self, aopdos_data, xlim=None, ylim=None):
        for name, data in aopdos_data.items():
            energy = data["Energy"]
            pdos = data["PDOS"]
            color = data["Color"]
    
            if pdos.shape[1] == 1:
                plt.fill_between(energy, np.ravel(pdos), color=color, alpha=0.5, label=name)
            else:
                pdosup = pdos[:, 0]
                pdosdw = -1 * pdos[:, 1]
                plt.fill_between(energy, pdosup, color=color, alpha=0.5, label=f"{name}")
                plt.fill_between(energy, pdosdw, color=color, alpha=0.5)
    
        plt.ylabel("PDOS (a.u.)", fontweight="bold")
        plt.xlabel("E - E$_f$ (eV)", fontweight="bold")
        plt.axvline(x=0, color='gray', linestyle='--')
        plt.yticks([])
        plt.xticks([])
        plt.legend()
    
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
    
        plt.show()

# Example usage:
# bs = BandEnergy()
# bs.get_bandstructure(data, fermi_energy, high_symmetry_points)
# bs.plot_shifted_band_structure(2, 0.5)
# 
# plotter = BandEnergy()
# files = {
#     "Br_s": ("blue", "../work/Pol_PDOS_Pure/(Br)_(s).dat"),
#     "Br_p": ("green", "../work/Pol_PDOS_Pure/(Br)_(p).dat"),
#     "Cs_s": ("orange", "../work/Pol_PDOS_Pure/(Cs)_(s).dat"),
#     "Cs_p": ("red", "../work/Pol_PDOS_Pure/(Cs)_(p).dat"),
#     "Cs_d": ("magenta", "../work/Pol_PDOS_Pure/(Cs)_(d).dat"),
#     "Pb_s": ("purple", "../work/Pol_PDOS_Pure/(Pb)_(s).dat"),
#     "Pb_p": ("brown", "../work/Pol_PDOS_Pure/(Pb)_(p).dat"),
#     "Pb_d": ("black", "../work/Pol_PDOS_Pure/(Pb)_(d).dat")
# }
# fermi_energy = 3.5926
# aopdos_data = plotter.get_aopdos(files, fermi_energy)
# plotter.plot_aopdos(aopdos_data, xlim=(-5, 5), ylim=(-80, 80))


class CombinedPlot:
    def __init__(self, band_energy):
        self.band_energy = band_energy

    def plot_combined(self, band_data, fermi_energy, kpoints, pdos_files, pdos_fermi_energy, band_xlim=None, band_ylim=None, pdos_xlim=None, pdos_ylim=None):
        # Create a figure with two subplots: one for band structure, one for PDOS
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), gridspec_kw={'width_ratios': [3, 2]})

        # Plot the band structure on the left subplot
        self.band_energy.data = band_data
        self.band_energy.fermi_energy = fermi_energy
        self.band_energy.k = np.unique(band_data[:, 0])
        self.band_energy.bands = np.reshape(band_data[:, 1], (-1, len(self.band_energy.k)))
        self.band_energy.high_symmetry_points = kpoints
        for i, band in enumerate(self.band_energy.bands):
            ax1.plot(self.band_energy.k, band - self.band_energy.fermi_energy, linewidth=1, alpha=0.5, color='r')
        ax1.set_xlim(min(self.band_energy.k), max(self.band_energy.k))
        if band_xlim is not None:
            ax1.set_xlim(band_xlim)
        ax1.set_ylim(-5, 5)
        if band_ylim is not None:
            ax1.set_ylim(band_ylim)
        if self.band_energy.high_symmetry_points is not None:
            x_coords = [x_coord for _, x_coord in self.band_energy.high_symmetry_points]
            x_labels = [symmetry_point for symmetry_point, _ in self.band_energy.high_symmetry_points]
            ax1.set_xticks(x_coords)
            ax1.set_xticklabels(x_labels)
            for _, x_coordinate in self.band_energy.high_symmetry_points:
                ax1.axvline(x=x_coordinate, color='k', linestyle='-', linewidth=1)
        ax1.set_ylabel("E-E$_f$ (eV)", fontweight='bold')
#        ax1.axvline(x=0, color='gray', linestyle='--')

        # Plot the PDOS on the right subplot
        aopdos_data = self.band_energy.get_aopdos(pdos_files, pdos_fermi_energy)
        for name, data in aopdos_data.items():
            energy = data["Energy"]
            pdos = data["PDOS"]
            color = data["Color"]
            if pdos.shape[1] == 1:
                ax2.fill_between(np.ravel(pdos), energy, color=color, alpha=0.5, label=name)
            else:
                pdosup = pdos[:, 0]
                pdosdw = -1 * pdos[:, 1]
                ax2.fill_between(energy, pdosup, color=color, alpha=0.5, label=f"{name}")
                ax2.fill_between(energy, pdosdw, color=color, alpha=0.5)
#        ax2.set_ylabel("PDOS (a.u.)", fontweight="bold")
        ax2.set_xlabel("PDOS (a.u.)", fontweight="bold")
#        ax2.axvline(x=0, color='gray', linestyle='--')
        ax2.set_yticks([])
        ax2.legend()
        if pdos_xlim is not None:
            ax2.set_xlim(pdos_xlim)
        if pdos_ylim is not None:
            ax2.set_ylim(pdos_ylim)

        plt.tight_layout()
        plt.show()

# Example usage:

# Assuming you have the band structure data, fermi energy, high symmetry points, and PDOS files defined

#band_data = np.array([[0.0, -1.0], [0.5, -0.5], [1.0, 0.0], [1.5, 0.5], [2.0, 1.0]])
#fermi_energy = 0.0
#high_symmetry_points = [("G", 0), ("X", 1), ("M", 2)]

# PDOS files dictionary
#files = {
#    "Br_s": ("blue", "../work/Pol_PDOS_Pure/(Br)_(s).dat"),
#    "Br_p": ("green", "../work/Pol_PDOS_Pure/(Br)_(p).dat"),
#    "Cs_s": ("orange", "../work/Pol_PDOS_Pure/(Cs)_(s).dat"),
#    "Cs_p": ("red", "../work/Pol_PDOS_Pure/(Cs)_(p).dat"),
#    "Cs_d": ("magenta", "../work/Pol_PDOS_Pure/(Cs)_(d).dat"),
#    "Pb_s": ("purple", "../work/Pol_PDOS_Pure/(Pb)_(s).dat"),
#    "Pb_p": ("brown", "../work/Pol_PDOS_Pure/(Pb)_(p).dat"),
#    "Pb_d": ("black", "../work/Pol_PDOS_Pure/(Pb)_(d).dat")
#}
#pdos_fermi_energy = 3.5926

# Create instances
#band_energy = BandEnergy()
#combined_plot = CombinedPlot(band_energy)

# Plot the combined figure
#combined_plot.plot_combined(band_data, fermi_energy, high_symmetry_points, files, pdos_fermi_energy, band_xlim=(0, 2), band_ylim=(-5, 5), pdos_xlim=(-5, 5), pdos_ylim=(-80, 80))
#

