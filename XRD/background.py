import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

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
