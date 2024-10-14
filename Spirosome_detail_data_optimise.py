# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:44:20 2024

@author: s2279999
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:29:13 2024

@author: s2279999
"""

from skimage.io import imread
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage import filters,measure
from skimage.filters import threshold_local








# Folders to analyse:
    
def load_image(toload):
    
    image=imread(toload)
    
    return image

def z_project(image):
    
    mean_int=np.mean(image,axis=0)
  
    return mean_int
def intensity_his_green(image,path,colour):
    pixel_values = image.flatten()
    plt.figure(figsize=(8,6))
    plt.hist(pixel_values, bins = 100, range =[0, 65535], color='green')
    plt.yscale('log')
    plt.xlabel('Pixel Intensity')
    plt.ylabel('Frequency')
    plt.title('Histogram of pixel intensity')
    plt.savefig(path+colour+'_histgram.png')
    plt.close()

def intensity_his_red(image,path,colour):
    pixel_values = image.flatten()
    plt.figure(figsize=(8,6))
    plt.hist(pixel_values, bins = 100, range = [0, 65535], color = 'red')
    plt.yscale('log')
    plt.xlabel('Pixel Intensity')
    plt.ylabel('Frequency')
    plt.title('Histogram of pixel intensity')
    plt.savefig(path+colour+'_histgram.png')
    plt.close()
    

# Subtract background:
def subtract_bg(image):
    background = threshold_local(image, 11, offset=np.percentile(image, 1), method='median')
    bg_corrected =image - background
    return bg_corrected

def threshold_image_std(input_image,thresh):
    threshold_value=filters.threshold_otsu(input_image)  
    # threshold_value=input_image.mean()+5*input_image.std()
    # threshold_value=thresh
    print(threshold_value)
    binary_image=input_image>threshold_value

    return threshold_value,binary_image,threshold_value

def threshold_image_standard(input_image,thresh):
     
    binary_image=input_image>thresh

    return binary_image

# Threshold image using otsu method and output the filtered image along with the threshold value applied:
    
def threshold_image_fixed(input_image,threshold_number):
    threshold_value=threshold_number   
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Label and count the features in the thresholded image:
def label_image(input_image):
    labelled_image=measure.label(input_image)
    number_of_features=labelled_image.max()
 
    return number_of_features,labelled_image
    
# Function to show the particular image:
def show(input_image,color=''):
    if(color=='Red'):
        plt.imshow(input_image,cmap="Reds")
        plt.show()
    elif(color=='Blue'):
        plt.imshow(input_image,cmap="Blues")
        plt.show()
    elif(color=='Green'):
        plt.imshow(input_image,cmap="Greens")
        plt.show()
    else:
        plt.imshow(input_image)
        plt.show() 
    
        
# Take a labelled image and the original image and measure intensities, sizes etc.
def analyse_labelled_image(labelled_image,original_image):
    measure_image=measure.regionprops_table(labelled_image,intensity_image=original_image,properties=('area','perimeter','centroid','orientation','major_axis_length','minor_axis_length','mean_intensity','max_intensity'))
    measure_dataframe=pd.DataFrame.from_dict(measure_image)
    return measure_dataframe

# This is to look at coincidence purely in terms of pixels

def coincidence_analysis_pixels(binary_image1,binary_image2):
    pixel_overlap_image=binary_image1&binary_image2         
    pixel_overlap_count=pixel_overlap_image.sum()
    pixel_fraction=pixel_overlap_image.sum()/binary_image1.sum()
    
    return pixel_overlap_image,pixel_overlap_count,pixel_fraction

# Look at coincidence in terms of features. Needs binary image input 

def feature_coincidence(binary_image1,binary_image2):
    number_of_features,labelled_image1=label_image(binary_image1)          # Labelled image is required for this analysis
    coincident_image=binary_image1 & binary_image2        # Find pixel overlap between the two images
    coincident_labels=labelled_image1*coincident_image   # This gives a coincident image with the pixels being equal to label
    coinc_list, coinc_pixels = np.unique(coincident_labels, return_counts=True)     # This counts number of unique occureences in the image
    # Now for some statistics
    total_labels=labelled_image1.max()
    total_labels_coinc=len(coinc_list)
    fraction_coinc=total_labels_coinc/total_labels
    
    # Now look at the fraction of overlap in each feature
    # First of all, count the number of unique occurances in original image
    label_list, label_pixels = np.unique(labelled_image1, return_counts=True)
    fract_pixels_overlap=[]
    for i in range(len(coinc_list)):
        overlap_pixels=coinc_pixels[i]
        label=coinc_list[i]
        total_pixels=label_pixels[label]
        fract=1.0*overlap_pixels/total_pixels
        fract_pixels_overlap.append(fract)
    
    
    # Generate the images
    coinc_list[0]=1000000   # First value is zero- don't want to count these. 
    coincident_features_image=np.isin(labelled_image1,coinc_list)   # Generates binary image only from labels in coinc list
    coinc_list[0]=0
    non_coincident_features_image=~np.isin(labelled_image1,coinc_list)  # Generates image only from numbers not in coinc list.
    
    return coinc_list,coinc_pixels,fraction_coinc,coincident_features_image,non_coincident_features_image,fract_pixels_overlap

# Rotate the image for chance
def rotate(matrix):
    temp_matrix = []
    column = len(matrix)-1
    for column in range(len(matrix)):
       temp = []
       for row in range(len(matrix)-1,-1,-1):
          temp.append(matrix[row][column])
       temp_matrix.append(temp)
    for i in range(len(matrix)):
       for j in range(len(matrix)):
          matrix[i][j] = temp_matrix[i][j]
    return matrix      

root_path="/Volumes/Tianxiao/20240826_gdl-c18 lipo_2024-08-26_16-40-30/"
# dirs_list = []
# for root,dirs,files in os.walk(root_path):
#     # if dirs != '':
#     if root =='/Volumes/Tianxiao/20240826_gdl-c18 lipo_2024-08-26_16-40-30/':
#         dirs_list.append(dirs)





# for xxxx in dirs_list:
#     for xxx in xxxx:







# These are the names of the files to image:

green_image="405_0.tif"
red_image="638_0.tif"

# Paths to analyse:

pathlist=[]

for i in range(1,4):
    for k in range(1,4):
        pathlist.append(root_path+"X0Y0R"+str(i)+"W"+str(k)+"_")

Output_all = pd.DataFrame(columns=['Number green','Number red','Number coincident','Number chance','Q','Green threshold','Red threshold'])
population_red_list = []
population_green_list = []
red_df_list = []
green_df_list = []
yellow_df_list = []
for path in pathlist:
    print('Current path is {}'.format(path))
    try:
        
      # Load the images
        green_path = path + green_image
        red_path = path + red_image
        green=load_image(green_path)
        red=load_image(red_path)
    
      # z-project - get the average intensity over the range. 
        
        green_flat=np.max(green,axis=0)
        
        red_flat=np.max(red,axis=0)
    
      # The excitation is not homogenous, and so need to subtract the background:
        
        green_bg_remove=subtract_bg(green_flat)
        population_green_list.append(green_bg_remove)
        intensity_his_green(green_bg_remove,path,green_image)
    
        red_bg_remove=subtract_bg(red_flat)
        population_red_list.append(red_bg_remove)
        intensity_his_red(red_bg_remove,path,red_image)
    
      # Threshold each channel: 
        
        thr_gr,green_binary,green_threshold_value = threshold_image_std(green_bg_remove,2844)
        
        thr_red,red_binary,red_threshold_value = threshold_image_std(red_bg_remove,3825)
       
      # Save the images 
        
        imsr = Image.fromarray(green_bg_remove)
        imsr.save(path+green_image+'_BG_Removed.tif')
        
        imsr = Image.fromarray(red_bg_remove)
        imsr.save(path+red_image+'_BG_Removed.tif')
        
       
        
      
        
        imsr = Image.fromarray(green_binary)
        imsr.save(path+green_image+'_Binary.tif')
        
        imsr = Image.fromarray(red_binary)
        imsr.save(path+red_image+'_Binary.tif')
        
      # Perform analysis 
       
        number_green,labelled_green=label_image(green_binary)
        print("%d feautres were detected in the green image."%number_green)
        measurements_green=analyse_labelled_image(labelled_green,green_flat)
        green_df_list.append(measurements_green)
           
        number_red,labelled_red=label_image(red_binary)
        print("%d feautres were detected in the red image."%number_red)
        measurements_red=analyse_labelled_image(labelled_red,red_flat)
        red_df_list.append(measurements_red)
      # Perform coincidence analysis
        
        green_coinc_list,green_coinc_pixels,green_fraction_coinc,green_coincident_features_image,green_non_coincident_features_image,green_fract_pixels_overlap=feature_coincidence(green_binary,red_binary)
        red_coinc_list,red_coinc_pixels,red_fraction_coinc,red_coincident_features_image,red_non_coincident_features_image,red_fract_pixels_overlap=feature_coincidence(red_binary,green_binary)
       
        number_of_coinc=len(green_coinc_list)
        
      # Need to account for chance due to high density
    
        green_binary_rot=rotate(green_binary) 
        
        chance_coinc_list,chance_coinc_pixels,chance_fraction_coinc,chance_coincident_features_image,chance_non_coincident_features_image,chance_fract_pixels_overlap=feature_coincidence(green_binary_rot,red_binary)
        
        number_of_chance=len(chance_coinc_list)
        
    # generate yellow coincident map
        
        yellow_map = green_coincident_features_image.astype(int)*green_flat/2 + red_coincident_features_image.astype(int)*red_flat 
        binary_yellow=threshold_image_standard(yellow_map,1)
        number_yellow, labelled_yellow = label_image(binary_yellow)
        measurements_yellow = analyse_labelled_image(labelled_yellow, yellow_map)
        yellow_df_list.append(measurements_yellow)
        
    #  Calculate an association quotient 
    
        Q=(number_of_coinc-number_of_chance)/(number_green+number_red-(number_of_coinc-number_of_chance))
        
    
        imsr = Image.fromarray(green_coincident_features_image)
        imsr.save(path+green_image+'_Coincident.tif')
       
        imsr = Image.fromarray(red_coincident_features_image)
        imsr.save(path+red_image+'_Coincident.tif')
        
        
        
    # Output
    
        Output_all = pd.concat([Output_all, pd.DataFrame({'Number green': [number_green],
                                                       'Number red': [number_red],
                                                       'Number coincident': [number_of_coinc],
                                                       'Number chance': [number_of_chance],
                                                       'Q': [Q],
                                                       'Green threshold':[green_threshold_value],
                                                       'Red threshold':[red_threshold_value]})], ignore_index=True)
    
        Output_all.to_csv(root_path+ '/All.csv', sep = '\t')
    except FileNotFoundError:
        continue
    except ValueError:
        continue
# write green, red and yellow events to csvs
combined_yellow = pd.concat(yellow_df_list, ignore_index = True)
combined_red = pd.concat(red_df_list, ignore_index = True)
combined_green = pd.concat(green_df_list, ignore_index = True)
combined_yellow.to_csv(root_path+'/Yellow_events.csv', index=True)
combined_green.to_csv(root_path+'/Green_events.csv', index=True)
combined_red.to_csv(root_path+'/Red_events.csv', index=True)

total_green = []
for arr in population_green_list:
    for element in arr:
        for data in element:
            total_green.append(data)
            

plt.figure(figsize=(8,6))
plt.hist(total_green, bins=100, range =[0, 65535],color='green')
plt.yscale('log')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.title('Histogram of pixel intensity')
plt.savefig(root_path+'/total_green_histgram.png')
plt.close()

total_red = []
for arr in population_red_list:
    for element in arr:
        for data in element:
            total_red.append(data)


plt.figure(figsize=(8,6))
plt.hist(total_red, bins=100, range =[0, 65535],color='red')
plt.yscale('log')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.title('Histogram of pixel intensity')
plt.savefig(root_path+'/total_red_histgram.png')
plt.close()
