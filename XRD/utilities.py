import numpy as np
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import os

class BackSub:
    def __init__(self, x, y):
        self.x = np.array(x)
        self.y = savgol_filter(np.array(y), 11, 3)
    
    def backsub(self, tol=1, show=False):
        L = len(self.y)
        lmda = int(0.50 * L / (self.x[0] - self.x[L-1]))  # Approx. # points for half width of peaks

        backsub_y = np.zeros(L)
        for i in range(L):
            if self.y[(i + lmda) % L] > tol * self.y[i]:  # Apply tolerance 'tol'
                backsub_y[(i + lmda) % L] = self.y[(i + lmda) % L] - self.y[i]
            else:
                if self.y[(i + lmda) % L] < self.y[i]:
                    backsub_y[(i + lmda) % L] = 0
        
        # Update y data with background-subtracted values
        self.y = backsub_y
        return self.x, backsub_y

# Example usage:
# x_data = [some x values]
# y_data = [some y values]
# bs = BackSub(x_data, y_data)
# x_sub, y_sub = bs.backsub(tol=1, show=True)


class peaks:
    def __init__(self, x, y):
        self.x = np.array(x)
        self.y = savgol_filter(np.array(y), 11, 3)
        self.peaks_array = None

    def find_peaks(self, height=None, distance=None, prominence=None, width=None):
        peaks, properties = find_peaks(self.y, height=height, distance=distance, prominence=prominence, width=width)
        
        # Store peaks as an array with each row [x_position, y_value]
        self.peaks_array = np.column_stack((self.x[peaks], self.y[peaks]))
        
        return self.peaks_array

# Example usage:
# x_data = np.linspace(0, 10, 100)
# y_data = np.sin(x_data) + 0.5 * np.random.randn(100)
# pf = peaks(x_data, y_data)
# peaks_array = pf.find_peaks(height=0.5)
# print(peaks_array)




class PlotTxtFiles:
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.data_files = self._get_txt_files()

    def _get_txt_files(self):
        """Get all .txt files in the specified folder."""
        return [f for f in os.listdir(self.folder_path) if f.endswith('.txt')]

    def read_file(self, file_name):
        """Read a single .txt file and return the data."""
        file_path = os.path.join(self.folder_path, file_name)
        try:
            data = np.genfromtxt(file_path, delimiter=',', skip_header=145, names=True)
        except ValueError as e:
            print(f"Error reading {file_name}: {e}")
            return None
        return data

    def plot_files(self):
        """Plot each .txt file in a separate subplot in an n x 1 matrix."""
        n = len(self.data_files)
        files_with_missing_data = []

        fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), sharex=True)

        # Ensure axes is always a list, even if n is 1
        if n == 1:
            axes = [axes]

        for i, file_name in enumerate(self.data_files):
            data = self.read_file(file_name)
            if data is None or 'Angle' not in data.dtype.names or 'PSD' not in data.dtype.names:
                files_with_missing_data.append(file_name)
                continue

            axes[i].plot(data['Angle'], data['PSD'], label=file_name[:-4])
            axes[i].set_ylabel('PSD')
            axes[i].legend()
            axes[i].grid(True)

        if files_with_missing_data:
            print("Files missing 'Angle' or 'PSD' columns:")
            for file_name in files_with_missing_data:
                print(file_name)

        axes[-1].set_xlabel('Angle')
        plt.suptitle('PSD vs Angle for Multiple Files')
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make space for the title
        plt.show()

# Example usage:
#folder_path = '../work/XRD/update_August_24/'
#plotter = PlotTxtFiles(folder_path)
#plotter.plot_files()
