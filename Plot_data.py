#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:02:20 2024

@author: tempchiao
"""

import os
import numpy as np
import random
import matplotlib.pyplot as plt
import csv
import statistics

# Define the working directory and file naming details
path = '/Volumes/Tianxiao/20240902_EVs/Fix/'
root = 'T30'

# Define parameters
file_numbers = 6               # Number of files collected
frequency = 10000              # Acquisition rate (in Hz)
acquisition_time = 30          # Duration of each file in seconds
display_range = 5              # Duration (in seconds) to visualize the data
green_threshold = 7            # Green channel threshold
red_threshold = 19             # Red channel threshold
sum_threshold = 20             # Sum threshold for event detection
threshold_type = 'and'         # Type of threshold ('sum' or 'and')

# Generate the list of file paths
file_list = []
for i in range(file_numbers):
    file_name = f"{root}_{i+1:02d}" if i > 0 else root
    file_list.append(os.path.join(path, file_name))

# Function to read data from a given file path
def read_data(file_path):
    data = np.genfromtxt(file_path, delimiter='\t')
    return data

# Function to generate a random range for data visualization
def generate_random_range():
    total_events = file_numbers * acquisition_time * frequency
    display_events = display_range * frequency
    start = random.randint(0, total_events - display_events)
    end = start + display_events
    return start, end

# Read data from all files and extract red and green channels
red_data = []
green_data = []
for file_path in file_list:
    data = read_data(file_path)
    for row in data:
        green_data.append(row[0])  # Assuming the first column is green
        red_data.append(row[-1])   # Assuming the last column is red

# Initialize counts and lists for events
red_count = 0
green_count = 0
FRET_count = 0
real_green = []
real_red = []

# Determine which thresholding method to use
if threshold_type == 'sum':
    for red, green in zip(red_data, green_data):
        if red + green > sum_threshold:
            FRET_count += 1
elif threshold_type == 'and':
    for red, green in zip(red_data, green_data):
        if red > red_threshold:
            real_red.append(red)
            red_count += 1
            if green > green_threshold:
                real_green.append(green)
                FRET_count += 1
        if green > green_threshold:
            green_count += 1

# Print results based on the threshold type
if threshold_type == 'sum':
    print(f'The number of FRET events is {FRET_count}')
else:
    print(f'The number of red events is {red_count}')
    print(f'The number of green events is {green_count}')
    print(f'The number of FRET events is {FRET_count}')

# Generate a random range for display
start, end = generate_random_range()

# Plotting the selected data range
negative_green_list = [-x for x in green_data[start:end]]
plt.figure(figsize=(15, 6))
plt.xlabel(f'{display_range}s')
plt.plot(range(start, end), negative_green_list, color='green', label='Green')
plt.plot(range(start, end), red_data[start:end], color='red', label='Red')
plt.tick_params(axis='y', which='major', labelsize=20)
plt.savefig(os.path.join(path, f'{root}_{display_range}s_preview.png'), transparent = True)

# Save event statistics to CSV
header = ['Red threshold', 'Green threshold', 'Red events', 'Green events', 
          'Red mean', 'Green mean', 'Red stdev', 'Green stdev', 'Red max', 'Green max']
body = [[red_threshold, green_threshold, red_count, green_count,
         statistics.mean(red_data), statistics.mean(green_data),
         statistics.stdev(red_data), statistics.stdev(green_data),
         max(red_data), max(green_data)]]
csv_file = os.path.join(path, f'{root}_number_of_events.csv')

# Append or create the CSV file
with open(csv_file, 'a', newline='') as file:
    writer = csv.writer(file)
    if not os.path.isfile(csv_file):  # If the file does not exist, write the header
        writer.writerow(header)
    writer.writerows(body)

# Display mean and standard deviation for Red and Green channels
print(f'Average intensity in Red channel: {np.mean(red_data):.2f}')
print(f'Standard deviation of intensity in Red channel: {np.std(red_data):.2f}')
print(f'Average intensity in Green channel: {np.mean(green_data):.2f}')
print(f'Standard deviation of intensity in Green channel: {np.std(green_data):.2f}')

# Uncomment below to calculate the average of real threshold events
# print('The average intensity of Green above Threshold is {:.3f}'.format(sum(real_green)/len(real_green)))
# print('The average intensity of Red above Threshold is {:.3f}'.format(sum(real_red)/len(real_red)))
