# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:23:46 2024

@author: s2279999
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import os

# The thresholds and filenames etc. need to be placed below

path = r'//Volumes/Tianxiao/20240410_EVs_CE/CEs/FRET/'


file_stem = "1_1_-10"  # This is the filename before the underscore.
csv_filename = path+'information.csv'
write = 'enabl'
grhis = 'enable'
rdhis = 'enable'
maxq = 'enabl'
number_of_files = 10

red_mean = 2.739694768
green_mean = 0.513702099
red_stdev = 9.473698325
green_stdev = 2.195933494

# Number of files in the folder (could make this automatic in the future).
channelA_thresh = 7    # Threshold for Channel A (Green).
channelB_thresh = 19  # Threshold for Channel B (Red).
channelA_AF = 0.6  # Autofluorescence
channelB_AF = 1.97
xtalk = 0.03  # Cross-talk from A to B

# Where channel A data will be stored


def load_files(number_of_files):
    channelA = []
    channelB = []  # Where channel B data will be stored
    for i in range(0, number_of_files):

        if(i < 1):
            filename = path + "/" + file_stem
        elif(i < 9):
            filename = path + "/" + file_stem + "_0" + str(i + 1)
        else:
            filename = path + "/" + file_stem + "_" + str(i + 1)
        a = 0  # Row counter
        with open(filename) as csvDataFile:  # Opens the file as a CSV
            csvReader = csv.reader(csvDataFile, delimiter='\t')  # Assigns the loaded CSV file to csvReader.
            for row in csvReader:
                channelA.append(row[0])  # For every row in in csvReader, the values are apended to green and red.
                channelB.append(row[1])
                a += 1

        print("Loaded %s, which contains %s rows." % (filename, a))

    channelA_arr = np.asarray(channelA, dtype=np.float32)  # Converts these to numpy arrays for vector calcs.
    channelB_arr = np.asarray(channelB, dtype=np.float32)
    return channelA_arr, channelB_arr

def histogram (colour, status):
    if status == 'enable':
        if colour == 'green':
            plt.figure(figsize=(8, 6))
            plt.hist(channelA_only_events,bins = 16, color ='green',rwidth=0.9, ec='black',label="Real Green Events")
            plt.yscale('log')
            plt.xlabel('Intensity')
            plt.ylabel('Number of count')
            plt.savefig(path + file_stem+ '_' +'Green_hist.pdf')
            
            return plt.show()
        if colour == 'red':
            plt.figure(figsize=(8, 6))
            plt.hist(channelB_only_events,bins = 16, color ='red',rwidth=0.9, ec='black',label="Real Red Events")
            plt.yscale('log')
            plt.xlabel('Intensity')
            plt.ylabel('Number of count')
            plt.savefig(path + file_stem+ '_' +'Red_hist.pdf')
            return plt.show()
    else:
        return 0

channelA_arr, channelB_arr = load_files(number_of_files)


# Now need to account for autofluorescence and crosstalk etc.

channelB_arr = (channelB_arr - xtalk * channelA_arr) - channelB_AF
channelA_arr = channelA_arr - channelA_AF

# This part is for the thresholding:

channelA_only_events = channelA_arr[(channelA_arr > channelA_thresh)]  # Total A events
channelB_only_events = channelB_arr[(channelB_arr > channelB_thresh)]  # Total B events
channelA_events = channelA_arr[np.logical_and(
    channelA_arr > channelA_thresh, channelB_arr > channelB_thresh)]  # A coincident events
channelB_events = channelB_arr[np.logical_and(
    channelA_arr > channelA_thresh, channelB_arr > channelB_thresh)]  # B coincident events

# Now need to account for chance events:

channelB_shuffle = channelB_arr.copy()
np.random.shuffle(channelB_shuffle)

channelA_chance = channelA_arr[np.logical_and(
    channelA_arr > channelA_thresh, channelB_shuffle > channelB_thresh)]  # These are the chance events
channelB_chance = channelB_shuffle[np.logical_and(
    channelA_arr > channelA_thresh, channelB_shuffle > channelB_thresh)]

# Now need to calculate Q value:

var_real_events = float(len(channelA_events))
var_A_events = float(len(channelA_only_events))
var_B_events = float(len(channelB_only_events))
var_chance_events = float(len(channelA_chance))
Q = float((var_real_events - var_chance_events) /
          (var_A_events + var_B_events - (var_real_events - var_chance_events)))

print('There were %s A events, %s B events, %s coincidence events, and %s chance events. Q = %f.' %
      (var_A_events, var_B_events, var_real_events, var_chance_events, Q))

# Now want histograms etc.

ln_events = np.log(channelB_events / channelA_events)
ln_chance = np.log(channelB_chance / channelA_chance)

textstr = 'Q = %.3f' % Q


plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = "12"
plt.figure(figsize=(8, 6))
plt.hist(ln_events, bins=60, range=[
         -3, 3], rwidth=0.9, ec='black', color='#ff0000', alpha=0.8, label="Real Events")
plt.hist(ln_chance, bins=60, range=[
         -3, 3], rwidth=0.9, ec='black', color='#cccccc', alpha=0.5, label="Chance Events")
plt.text(0.05, 0.90, textstr, transform=plt.gca().transAxes)
plt.legend(loc='upper right')
plt.xlabel('Z=$ln(I_{B}/I_{A}$)')
plt.ylabel('Number of events')
plt.savefig(path + file_stem + '_' + 'Mix_lnz.pdf')
plt.show()

FRET = channelB_events / (channelA_events + channelB_events)
plt.hist(FRET, bins=20, range=[
         0, 1], rwidth=0.9, ec='black', color='#ff0000', alpha=0.8, label="Real Events")
plt.text(0.05, 0.90, textstr, transform=plt.gca().transAxes)

FRET_ave = np.mean(FRET)
# print(FRET_ave)
plt.xlabel('Proximity Ratio')
plt.ylabel('Number of events')
plt.savefig(path + file_stem+ '_' +'lnz.pdf')
plt.show()


def maxQ(status):
    if status == 'enable':
        q_vals = np.zeros(shape=(20, 20))

        for A in range(20):
            for B in range(20):
                channelA_only_events = channelA_arr[(
                    channelA_arr > A)]  # Total A events
                channelB_only_events = channelB_arr[(
                    channelB_arr > B)]  # Total B events
                channelA_events = channelA_arr[np.logical_and(
                    channelA_arr > A, channelB_arr > B)]  # A coincident events

                # Now need to account for chance events:

                channelB_shuffle = channelB_arr.copy()
                np.random.shuffle(channelB_shuffle)

                channelA_chance = channelA_arr[np.logical_and(
                    channelA_arr > A, channelB_shuffle > B)]  # These are the chance events

                # Now need to calculate Q value:

                var_real_events = float(len(channelA_events))
                var_A_events = float(len(channelA_only_events))
                var_B_events = float(len(channelB_only_events))
                var_chance_events = float(len(channelA_chance))
                Q = float((var_real_events - var_chance_events) /
                        (var_A_events + var_B_events - (var_real_events - var_chance_events)))

                q_vals[A][B] = Q

        maximum_Q = np.amax(q_vals)
        result = np.where(q_vals == np.amax(q_vals))
        ThresholdA, ThresholdB = result

        print ('The maximum value of Q is %.3f, with a threshold of %s in channel A, and %s in channel B.'%(maximum_Q,str(ThresholdA),ThresholdB))

        contourplot = plt.contourf(q_vals, 20, origin='lower')
        cbar = plt.colorbar(contourplot)
        plt.xlabel("Channel B Threshold")
        plt.ylabel("Channel A Threshold")
        cbar.ax.set_ylabel('Q')
        global path
        global file_stem
        plt.savefig(path + file_stem+ '_' +'MaxQ.pdf')
        plt.show()
    else:
        return 0

if write == 'enable':
    header = ['Green threshold','Red thershold','Green events','Red events','Number of Coincidence','Number of chance','Q']
    
    file_exists = os.path.isfile(csv_filename)
    
    with open (csv_filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        if not file_exists:
            writer.writerow(header)
        
        row1 = [float(channelA_thresh),
                float(channelB_thresh),
                float(var_A_events),
                float(var_B_events),
                float(var_real_events),
                float(var_chance_events),
                float(Q)]
        row2 = [file_stem]
        writer.writerow(row2)
        writer.writerow(row1)


# plt.figure(figsize=(8, 6))
# plt.hist(channelA_only_events,bins = 16, range = (0,200),color ='green',rwidth=0.9, ec='black',label="Real Green Events")
# plt.yscale('log')
# plt.xlabel('Intensity')
# plt.ylabel('Number of count')
# plt.savefig(path + file_stem+ '_' +'Green_hist.pdf')
histogram('green',grhis)
histogram('red',rdhis)
""" plt.figure(figsize=(8, 6))
plt.hist(channelB_only_events,bins = 16, color ='red',rwidth=0.9, ec='black',label="Real Red Events")
plt.yscale('log')
plt.xlabel('Intensity')
plt.ylabel('Number of count')
plt.savefig(path + file_stem+ '_' +'Red_hist.pdf') """
maxQ(maxq)
