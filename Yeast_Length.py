#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:26:47 2018

@author: Rory Maizels
"""


### IMPORTS ###

"""
Program for segmentation and object analysis of S. cerevisiae microscopy images.
Specifically used to determine average length (in pixels) of cells in image,
used for cell-cycle research.

Lets user choose .tif or containing the set of images to be analysed. Uses
canny edge detection to detect cell objects, excludes objects that have
unrealistic size or shape (determined by size, circularity and convexity).
"""

import numpy as np
import scipy.stats as st
import skimage
from skimage.viewer import ImageViewer
from skimage.feature import canny
import os
from skimage import io
from skimage.exposure import rescale_intensity
from skimage.measure import label, regionprops
import matplotlib.pyplot as plt
from skimage import morphology
from scipy import ndimage as ndi
from skimage.color import label2rgb
from math import pi
from tkinter import filedialog
import tkinter as tk
import sys
import time

### END OF MODULES ###

### FUNCTIONS ###

def load_file():
    """
    let user select tif file.
    """
    root = tk.Tk()
    filename = filedialog.askopenfilename()
    root.destroy()
    return filename

def analyse_images(images):
    """
    given images, segment each image and calculate lengths of cell objects
    return list of all lengths, list of lengths per image, mean length per image,
    standard deviation per image, segmentation images, total average, total
    standard deviation and 95% confidence interval for this.
    """
    
    lollen = []
    alllen = []
    means = []
    sds = []
    lolabov = []
    try:
        number = len(images[:,0,0])
    except IndexError:
        number = 1
    for i in range(number):
        image = images[i,:,:]
        edges = canny(image/255., sigma = 0.9, low_threshold= 0.8, high_threshold= 0.85,
                      use_quantiles=True)
        fill = ndi.binary_fill_holes(edges)
        nosmall = morphology.remove_small_objects(fill,2500)
        invfill = morphology.remove_small_objects(fill,10000)
        nobig = np.logical_xor(nosmall,invfill)
        labeled = label(nobig)
        imlabov = label2rgb(labeled, image=nobig)
        lolabov.append(imlabov)
        props = regionprops(labeled)
        length = []
        for prop in props:
            A = prop['area']
            P = prop['perimeter']
            circularity = 4*pi*(A/P**2)
            if circularity > 0.4:
                S = prop['solidity']
                if S > 0.9:
                    length.append(prop['major_axis_length'])
                    alllen.append(prop['major_axis_length'])
        mean = np.mean(length)
        sd = np.std(length)
        lollen.append(length)
        means.append(mean)
        sds.append(sd)
    totav = np.mean(alllen)
    totsd = np.std(alllen)
    confi = st.t.interval(0.95, len(alllen)-1, loc=np.mean(alllen),
                       scale=st.sem(alllen))
    return lollen, alllen, means, sds, lolabov, totav, totsd, confi

def spacer():
    """
    just to space out user dialog screen.
    """
    print()
    print()
    time.sleep(0.2)
    print("...")
    print()
    print()

        
### END OF FUNCTIONS ###       
        
### PROGRAM ###
             
filename = load_file()
images = io.imread(filename)
images = rescale_intensity(images) 

(list_of_lengths, all_lengths, image_means, 
 image_sds, seg_images, total_ave, total_std, ci) = analyse_images(images)
options = ['1','2','3','4','5','6']
choice = '0'
while choice != '6':
    print("""1. individual cell lengths segregated by image.
2. all individual cell lengths in one list.
3. mean and SD values of each image.
4. total average, SD and statistics of all cells measured.
5. visual plot of each image segmentation.
6. exit program.""")
    while choice not in options:
        choice = input("Select an Option: ")
        spacer()
    if choice == '1':
        for i in range(len(list_of_lengths)):
            print("Cell Image " + str(i+1) + ":")
            print(list_of_lengths[i])
            choice = '0'
            spacer()
    if choice == '2':
        print("All Cell Lengths: ")
        print(all_lengths)
        choice = '0'
        spacer()
    if choice == '3':
        print("Image Means: ")
        print(image_means)
        print()
        print("Image Standard Deviations: ")
        print(image_sds)
        choice = '0'
        spacer()
    if choice == '4':
        n = len(all_lengths)
        print("Total Mean: ")
        print(total_ave)
        print()
        print("Total Standard Deviation: ")
        print(total_std)
        print()
        print("95% Confidence Interval: ")
        print(ci)
        print()
        print("n =" + str(n))
        spacer()
        choice = '0'
    if choice == '5':
        num = len(list_of_lengths)
        plt.figure(figsize=(9, 4*num))
        for i in range(num):
            plt.subplot(num,2,(2*i+1))
            plt.imshow(images[i,:,:],cmap='gray')
            plt.title("Original (Image " + str(i+1) + ")")
            plt.subplot(num,2,(2*i+2))
            plt.imshow(seg_images[i])
            plt.title("Segmented (Image " + str(i+1) + ")")
        plt.show()
        spacer()
        print("""NOTE: Some objects in the segmented image will be excluded
from calculations.""")
        spacer()
        choice = '0'
    if choice == '6':    
        print("bye bye!")
    
    




