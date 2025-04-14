# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 21:04:42 2025

@author: ychuang
"""


from SOLPS_input.header import *


if __name__ == "__main__":


     # Sample data
     x = [1, 2, 3, 4, 5]
     y = [2, 4, 1, 8, 7]

     # Create the plot
     plt.figure()
     plt.plot(x, y, label="Line", marker='o', color='blue')

     # Add title and labels
     plt.title("Simple Line Plot")
     plt.xlabel("X Axis")
     plt.ylabel("Y Axis")

     # Add legend
     plt.legend()
     plt.savefig("test.png", format="png")

     # Show the plot
     plt.show()