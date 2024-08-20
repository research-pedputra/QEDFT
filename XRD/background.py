import numpy as np
import matplotlib.pyplot as plt

class BackSub:
    def __init__(self, x, y):
        """
        Initialize the BackSub object with x and y data.

        Parameters
        ----------
        x : array-like
            The x-axis data.
        y : array-like
            The y-axis data.
        """
        self.x = np.array(x)
        self.y = np.array(y)
    
    def backsub(self, tol=1, show=False):
        """
        Perform background subtraction operation.

        Parameters
        ----------
        tol : float, optional
            Tolerance for background subtraction.
            Default value is 1.
        show : bool, optional
            Whether to show the plot of the XRD chart.
            Default value is False.
        
        Returns
        -------
        tuple
            A tuple containing the x data and the background-subtracted y data.
        """
        
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

        if show:
            plt.plot(self.x, self.y)
            plt.show()

        return self.x, backsub_y

# Example usage:
# x_data = [some x values]
# y_data = [some y values]
# bs = BackSub(x_data, y_data)
# x_sub, y_sub = bs.backsub(tol=1, show=True)
