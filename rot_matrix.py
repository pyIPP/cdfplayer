import numpy as np


class ROT_MATRIX:


    def __init__(self, alpha, x_in, y_in, x_cen=0, y_cen=0):

        x_in2  = x_in - x_cen      
        y_in2  = y_in - y_cen
        x_out2 = x_in2*np.cos(alpha) - y_in2*np.sin(alpha)
        y_out2 = x_in2*np.sin(alpha) + y_in2*np.cos(alpha)
        self.x = x_out2 + x_cen
        self.y = y_out2 + y_cen
