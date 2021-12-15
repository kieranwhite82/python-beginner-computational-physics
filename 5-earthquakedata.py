# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:26:31 2021

@author: uqkwhi18

Program that requests the user to input a month and day, and responds by
reading in a datafile on earthquakes in 2014 with the correct format for each
record, processing it and printing out whether there was a major earthquake on
that day.

Input
Month and day.

Output
Date.
If there was a major earthquake.
What the magnitude of the earthquake/s was/were.
Longitude and latitude of earthquake.
The number of earthquakes that occured on the entered date.
"""
import numpy as np

# Define path to data file.
path = 'C:/Users/kiera/OneDrive/Documents/Job/BRC/Coding/worldearthquakes20'\
       '14.dat'

# Set the cloumn widths in preparation for reading in the data file: delimiters
# are variable.
# Column locations are corrected for text reader starting at 1 and python
# starting at 0.
column_locations = np.array([1, 6, 9, 15, 26, 38, 48, 52])
column_locations = column_locations - 1
widths = column_locations[1:] - column_locations[:-1]

# Read in data file defined by path variable.
data = np.genfromtxt(path, dtype=float, delimiter=widths, autostrip=True)

# Get user input on the date: day and month, in number form.
inputm = eval(input("Enter the month date number in which to check a day for \
major earthquakes (August = 8, December=12): "))
inputd = eval(input("Enter the day date number on which to check for major \
earthquakes: "))

# Print out date and column titles for data on earthquakes that occured on the
# user entered date.
print("\nThe date the following earthquakes occured on is {}/{}/2014.".format(\
      inputd, inputm))
print("\nMagnitude\tLatitude\tLongitude")

# Initialise integers for looping (i) and counting (j).
i = 0
j = 0

# Loop to read over each row in the data. If the day and month are both equal
# to the user-entered values, the data on the earthquakes is printed.
for i in range(0, len(data[:, 1]) - 1):
    if data[i, 1] == inputm and data[i, 2] == inputd:
        print("{:.2f}\t\t{:.2f}\t\t{:.2f}".format(data[i, 6], data[i, 3], data\
                                                  [i, 4]))
        j += 1
    else:
        i += 1

# Print the number of earthquakes that occured on the user-entered date.
print("\nThe number of earthquakes that occured on {}/{}/2014 was {}.".format(
       inputd, inputm, j))
