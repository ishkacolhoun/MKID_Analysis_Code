#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:21:41 2022

@author: ishka
"""
#This code plots vatious plots for a single resonator at various power levels.

#note, works for IQ files saved using fsweep macro saved on the VNA using the VNA code IQSweepArray_Final_version.vbs
#its important to use this to save the data as the naming format is used in this code to find the power of the resonators.


#import usual packages 
import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import argrelmin
from scipy import optimize
from scipy.optimize import curve_fit
import pylab
import matplotlib.pyplot as plt
import array as arr
import csv 
from matplotlib.backends.backend_pdf import PdfPages
from datetime import date 
import os

# Input the name of the resonator you want to plot
name_of_resonator ="3435960000"

# This is the multiples of FWHM you want to plot from the min freq, anything outside this bound is trimmed from the data for the purpose of plots
distmult = 4

#This is just used in the figure plots to make sure to plot different plots for the Raw and Fitted IQ DATA for a single power
n = 10

#creates an empty list for the resonator files
resonator_files = []

#creates a list with all the files for a certain resonator, the directory here should contain all .csv files you want to plot, note they can be in subdirectories
for path, currentDirectory, files in os.walk("/home/ishka/Desktop/Library/MKIDS/Current Resonator"):
    for file in files:
        if file.endswith(name_of_resonator +".csv"):
            resonator_files.append(os.path.join(path, file))


#gets todays date. This is used in the name of the results files that are saved. 
today = date.today() 
todaystr = today.strftime("%Y_%m_%d") 

#E.17 from Gao thesis - will be used later to  amplitude data fit curve
def amplitudeequation(f, A1, A2, A3, A4, Qr, fr):
    return A1 + A2*(f - fr) + ((A3 + A4*(f - fr)) / (1 + 4 * Qr**2  * ((f - fr)/fr)**2 ))  

#E.11 from Gao thesis - can be used to fit phase data
def phaseequation(f, theta0, Qr, fr):
    return -theta0 + 2*np.arctan(2*Qr*(1 - (f/fr))) 

#E.11 again just in degrees this time
def phaseequation_degrees(f, theta0, Qr, fr):
    return(-theta0 + 2*np.arctan(2*Qr*(1 - (f/fr)))) * 180/np.pi

#The directory the plot pdfs will be saved
path = "/home/ishka/Desktop/Library/MKIDS/Current Resonator/Plots"  



#Sets the initial limit in range for the frequency
xlim_right = 0
xlim_left = 1000000000000


#creates empty lists this is used to calculate the average Q values later on

lst_of_Q_phase = []
lst_of_Qi_phase = []
lst_of_Qc_phase = []
lst_of_fr_phase = []

#creates the list where the power of the resonators will be saved this is a reason the naming structure of the sweep files is important
power_list = [] 


#loops through plotting code for each resonator file

for fullname in resonator_files:
    Power = fullname.replace( "_" + name_of_resonator + ".csv","") # retrieves the power from the name of the resonator file
    power = - int(Power[-2] + Power[-1])
    power_list.append(power)
    power_list.sort()
    print("Power:",power)
    
    #start reader and calculate length of file
    with open(fullname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        row_count = sum(1 for row in csv_reader)  # fileObject is your csv.reader  
        
    number_of_values = row_count - 10 #accounts for extra lines at start and end of file
    
    #declare 3 arrays, for frequency and I and Q values
    freqs = np.zeros(number_of_values)      
    Ivalues = np.zeros(number_of_values)
    Qvalues = np.zeros(number_of_values)

    
    #open reader again, this time to read in data 
    with open(fullname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        
        line_count = 0
        for row in csv_reader:
            if line_count <= 6: #opening lines, before actual data 
                line_count += 1
    
            elif line_count < row_count - 1 : #this is the actual IQ data
                freqs[line_count-7] = row[0] #saves data to array from csv file
                Ivalues[line_count-7] = row[1] #saves data to array from csv file
                Qvalues[line_count-7] = row[2]  #saves data to array from csv file
                
                line_count += 1
            if line_count >= (number_of_values + 7):
                break

    freqsMHz = freqs / 10**6
    
    
    #calculate amplitude of IQ data 
    amplitudes = np.sqrt((Ivalues)**2 + (Qvalues)**2)
    
    #find phases of raw data relative to origin
    #have to do this for 4 different quadrants    
    phases = np.arctan2(Qvalues, Ivalues)
    
    
    normalizedamplitudes = amplitudes / np.max(amplitudes) #normalizes magnitude data relative to max value
    squaremagnitudes = normalizedamplitudes**2 #calculates square magnitude 
    
    p0 = np.array([1, 1, 1, 1, 10000000, np.mean(freqsMHz)]) #initial guesses for fit
    
    #need to define sigma array to weight certain points heavier 
    resonant_index = np.argmin(abs(squaremagnitudes)) #index of resonant point in array
    
    weight_amplitude = np.ones(number_of_values) #set relative error everywhere to 1
    weight_amplitude[resonant_index - 25 : resonant_index + 25] = 0.15 #set error around resonant to smaller 
    
    
    #fits square magnitude data to amplitude equation
    #popt is the array of fitted parameters
    #pcov is the covariance of popt
    popt, pcov = curve_fit(amplitudeequation, freqsMHz, squaremagnitudes, p0, sigma = weight_amplitude, absolute_sigma=False, maxfev=10000)   
    
    
    Q = popt[4] #total Q is taken from the fit
    Q = abs(Q)
    
    fr = popt[5]  #resonant frequency is taken from the fit
    
    Qi = Q / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 
    Qi = abs(Qi)
    
    Qc = Q / (1 - np.min(normalizedamplitudes)) #calulculates coupling Q using equation 29 from Zmuidzinas review
    Qc = abs(Qc)
    
    
    squaremagnitudesdb = 10*np.log10(squaremagnitudes) #converts square magnitude data to dbs
    
    #calculates  approx FWHM to use in trimming the edges of the data
    
    def closest(array, K): # function to return the index of the closest element to K in an array
      
     idx = (np.abs(array - K)).argmin()
     return idx
    # this code calculates a rough Full width Half Max value
    plot_array = squaremagnitudesdb
    half_min = min(plot_array)/2
    min_index = (np.ndarray.tolist(plot_array)).index(min(plot_array))
    
    # find index of element closest to the half max value
    closest_indx = np.abs(plot_array - half_min).argmin()
    
    # the distance from the min freq value that we will set as our left and right most limits for freq
    distance = 2 * distmult * np.abs(freqsMHz[closest_indx]-freqsMHz[min_index])
    
    #defining the closest left and right index
    closest_indx_left = closest(freqsMHz, freqsMHz[min_index] - distance)
    closest_indx_right = closest(freqsMHz, freqsMHz[min_index] + distance)
    
   
    # defines new freq I and Q values with edges trimmed around the left and right limit
    freqsMHz_new = freqsMHz[closest_indx_left:closest_indx_right]
    Ivalues_new = Ivalues[closest_indx_left:closest_indx_right]
    Qvalues_new = Qvalues[closest_indx_left:closest_indx_right]
    
    number_of_values_new = len(freqsMHz_new)
    
    
# had to redo lines 98 to 165 to make them use the new trimmed data we have just created
    





     
    #calculate amplitude of IQ data 
    amplitudes = np.sqrt((Ivalues_new)**2 + (Qvalues_new)**2)
    
    #find phases of raw data relative to origin
    #have to do this for 4 different quadrants    
    phases = np.arctan2(Qvalues_new, Ivalues_new)
    
    
    normalizedamplitudes = amplitudes / np.max(amplitudes) #normalizes magnitude data relative to max value
    squaremagnitudes = normalizedamplitudes**2 #calculates square magnitude 
    
    p0 = np.array([1, 1, 1, 1, 1000000, np.mean(freqsMHz_new)]) #initial guesses for fit
    
    #need to define sigma array to weight certain points heavier 
    resonant_index = np.argmin(abs(squaremagnitudes)) #index of resonant point in array
    
    weight_amplitude = np.ones(number_of_values_new) #set relative error everywhere to 1
    weight_amplitude[resonant_index - 25 : resonant_index + 25] = 0.15 #set error around resonant to smaller 
    
    
    #fits square magnitude data to amplitude equation
    #popt is the array of fitted parameters
    #pcov is the covariance of popt
    popt, pcov = curve_fit(amplitudeequation, freqsMHz_new, squaremagnitudes, p0, sigma = weight_amplitude, absolute_sigma=False, maxfev=10000)   
    
    
    Q = popt[4] #total Q is taken from the fit
    Q = abs(Q)
    
    fr = popt[5]  #resonant frequency is taken from the fit
    
    Qi = Q / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 
    Qi = abs(Qi)
    
    Qc = Q / (1 - np.min(normalizedamplitudes)) #calulculates coupling Q using equation 29 from Zmuidzinas review
    Qc = abs(Qc)
    
    n = n + 1
    
    
    squaremagnitudesdb = 10*np.log10(squaremagnitudes) #converts square magnitude data to dbs
    
    
    
    
    
    
    
    #next need to convert this to IQ data and plot versus raw data 
    squaremagnitudefit = amplitudeequation(freqsMHz_new, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) #saves square magnitude fit data
    
    fitamplitude = np.sqrt(np.abs(squaremagnitudefit)) #takes sqrt to get amplitude
    
    fitamplitudescaled = fitamplitude * np.max(amplitudes) #unnormalizes data relative to max value
     
    
    
    #get centre of loop
    #find centre of loop
    raw_data_xc = (np.max(Ivalues_new) + np.min(Ivalues_new)) / 2
    raw_data_yc = (np.max(Qvalues_new) + np.min(Qvalues_new)) / 2
    
    
    #Next step, remove cable delay term from data. 
    #To do this, normalize by cable delay loop
    #For now, using circle with radius 4, centered on origin
    #Will have to use actual data, after next cooldown
    cabledelayphase = np.zeros(number_of_values_new)
    
    #first find total angle, between first and last point 
    #find initial phase
    
    initial_phase = np.arctan2(Qvalues_new[0], Ivalues_new[0])
    final_phase = np.arctan2(Qvalues_new[number_of_values_new-1] , Ivalues_new[number_of_values_new-1])
    
    total_phase = abs(final_phase - initial_phase)
    
    phase_per_freq = total_phase / (number_of_values_new - 1)
    
    i = 0
    for x in cabledelayphase:    
        cabledelayphase[i] = initial_phase - phase_per_freq*i        
        i += 1
    
    
    #need to find amplitude of cable delay loop - take average of start and end  of loop
    amplitudes = np.sqrt(Ivalues_new**2 + Qvalues_new**2)
        
    cableamplitude = np.max(amplitudes) #takes max value 
    cabledelayI = cableamplitude*np.cos(cabledelayphase) #finds I component of cable
    cabledelayQ = cableamplitude*np.sin(cabledelayphase) #finds Q component of cable
        
    
    #substract cable delay from IQ data to normalize
    Inormalized = Ivalues_new - cabledelayI #chceck if this is correct?
    Qnormalized = Qvalues_new - cabledelayQ
    
    #amplitudes_normalized = np.zeros(number_of_values)
    
    amplitudes_normalized = np.sqrt(Inormalized**2 + Qnormalized**2)
    
    #next have to move loop to be centred on (0, 0)
    
    xc = (np.max(Inormalized) + np.min(Inormalized)) / 2 
    yc = (np.max(Qnormalized) + np.min(Qnormalized)) / 2
    
    #find alpha - the argument of zc, the centre of the normalized circle
    alpha = np.arctan2(yc, xc)
    
    translatedcircleI = Inormalized - xc #finds I component of cable
    translatedcircleQ = Qnormalized - yc #finds Q component of cable
        
    amplitudes_translated = np.sqrt(translatedcircleI**2 + translatedcircleQ**2)
    
    
    theta = np.arctan2(translatedcircleQ, translatedcircleI)
    r_translated = np.sqrt((translatedcircleI)**2 + (translatedcircleQ)**2)
    
    #now rotate everything by alpha 
    
    #need to rotate every point in this circle by -alpha 
    theta_rotated = theta - alpha #rotates by alpha
    I_rotated = r_translated*np.cos(theta_rotated) #finds I component of cable
    Q_rotated = r_translated*np.sin(theta_rotated) #finds Q component of cable
        
    #amplitudes_rotated = np.zeros(number_of_values)
    amplitudes_rotated = np.sqrt(I_rotated**2 + Q_rotated**2)
    
    thetafinal = np.unwrap(theta_rotated)
    
    #the following shifts theta so that it is centered on 0 deg
    
    midpt_thetafinal = (max(thetafinal) - min(thetafinal)) / 2 + min(thetafinal) #finds the offset needed for the midpoint to be at 0 deg
    for n1 in range(len(thetafinal)): # shifts every element in thetafinal by the offset midpt_thetafinal 
        thetafinal[n1] = thetafinal[n1] - midpt_thetafinal
    
    
    #now, fit this to phase equation 
    p_phase = np.array([1, Q, fr])
      
    popt_phase, pcov_phase = curve_fit(phaseequation, freqsMHz_new, thetafinal, p_phase)

    
    theta0_phase = popt_phase[0]
    Qr_phase = popt_phase[1]
    fr_phase = popt_phase[2]

    
    Q_phase = abs(Qr_phase)
    
    Qi_phase = Q_phase / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 

    Qi_phase = abs(Qi_phase)
 
    
    Qc_phase = Q_phase / (1 - np.min(normalizedamplitudes)) #calulculates instrinsic Q using equation 29 from Zmuidzinas review
    Qc_phase = abs(Qc_phase)
    
    
    #final step, convert phase and amplitude fit back to IQ data and replot fit
    fitted_phase = phaseequation(freqsMHz_new, popt_phase[0], popt_phase[1], popt_phase[2])
    
    #first step, get I and Q data and rotate and translate back to position 
    
    
    #fitted_phase_rotated = np.zeros(number_of_values)
    fitted_phase_rotated = fitted_phase + alpha
    
    fitI = np.mean(amplitudes_rotated)*np.cos(fitted_phase) #finds I component of fit
    fitQ = np.mean(amplitudes_rotated)*np.sin(fitted_phase) #finds Q component of fit
            
    fitIrotated = np.mean(amplitudes_rotated)*np.cos(fitted_phase_rotated) #finds I component of fit
    fitQrotated = np.mean(amplitudes_rotated)*np.sin(fitted_phase_rotated) #finds Q component of fit
        
    fitItranslated = fitIrotated + xc
    fitQtranslated = fitQrotated + yc
        
    fitIunnormalized = fitItranslated + cabledelayI
    fitQunnormalized = fitQtranslated + cabledelayQ    
    
    #get phase of unnormalize data
    fitphaseunnormalized = np.arctan2(fitQunnormalized, fitIunnormalized)    
        
    #last thing to do, combine fitted amplitude and phase data    
    Ifitfinal = fitamplitudescaled*np.cos(fitphaseunnormalized) #finds I component of cable
    Qfitfinal = fitamplitudescaled*np.sin(fitphaseunnormalized) #finds Q component of cable
    
    #functions to calculate the average Q values
    def Average_Q_phase():
        return sum(lst_of_Q_phase) / len(lst_of_Q_phase)
    lst_of_Q_phase.append(Q_phase)
    
    def Average_Qi_phase():
        return sum(lst_of_Qi_phase) / len(lst_of_Qi_phase)
    lst_of_Qi_phase.append(Qi_phase)
    
    def Average_Qc_phase():
        return sum(lst_of_Qc_phase) / len(lst_of_Qc_phase)
    lst_of_Qc_phase.append(Qc_phase)
    
    def Average_fr_phase():
        return sum(lst_of_fr_phase) / len(lst_of_fr_phase)
    lst_of_fr_phase.append(fr_phase)
    
    # the below calculates the value of r^2 of the fits
        # residual sum of squares
    ss_res = np.sum((Qvalues_new - Qfitfinal) ** 2)
    
    # total sum of squares
    ss_tot = np.sum((Qvalues_new - np.mean(Qvalues_new)) ** 2)
    
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    r2_formatted = "{:.4f}".format(r2)
    print('r^2 of fit = ' + r2_formatted)
    
    #the name of the saved pdf
    savefilename = path + "/" + name_of_resonator + "_Fit_" + todaystr + '.pdf'

    
    #creates pdf to save plots
    pp = PdfPages(savefilename)       #commentig this out  while fixing to fit phase too
  
    
  
    #now plot just origianl IQ data and fitted data for all powers on one graph
 
    plt.figure(1)
    plt.scatter(Ivalues_new, Qvalues_new, marker='.', s=0.01)
    plt.plot(Ifitfinal, Qfitfinal, linewidth = 1)
    plt.xlabel('I')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Q')
    plt.title('Raw and Fitted IQ Data')
    plt.grid()
    
    plt.savefig(pp, format='pdf')
    
 # plots the raw and fitted IQ for each resonator
    
    plt.figure(n)
    plt.scatter(Ivalues_new, Qvalues_new, marker='.', s=0.01)
    plt.plot(Ifitfinal, Qfitfinal, linewidth = 1)
    plt.xlabel('I')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Q')
    plt.title('Raw and Fitted IQ Data resonator @ '+ str(power) + 'db')
    plt.grid()   
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
    boxstr = f'Avrg_Q = {round(int(Average_Q_phase()), -2)} \n Avrg_Qc = {round(int(Average_Qc_phase()), -2)} \n Avrg_Qi = {round(int(Average_Qi_phase()), -2)} \n Avrg_fr = {round(Average_fr_phase(), 4)} MHz \n R squared = {r2_formatted}' #contents of text box
    plt.text(min(Ivalues_new) ,min(Qvalues_new),boxstr, fontsize=11, bbox=boxprops)#adds textbox to plot
      
    plt.savefig(pp, format='pdf')
    
    
    #plots just origianl IQ but only every 10th point 
    plt.figure(2)
    nth_Ivalue = Ivalues_new[::10]
    nth_Qvalue = Qvalues_new[::10]
    plt.scatter(nth_Ivalue, nth_Qvalue, marker='.', s=0.05)
    plt.xlabel('I')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Q')
    plt.title('Raw IQ Data(every 10th point) ')
    plt.grid()
  
    plt.savefig(pp, format='pdf')
    

    #plots normalised S21 vs Freq for both raw and fitted data and all powers on the same graph    

    plot_label = "Resonator Power: " + str(power) + "db"
    #plots square magnitudes verus frequency in dbs 
    plt.figure(3)
    plt.scatter(freqsMHz_new, squaremagnitudesdb, marker='.', s=0.005, label = plot_label)
    plt.title('Normalized $|S_{21}|^2$ (dB) vs Frequency of Raw Data and Fit')
    plt.legend(loc=2, prop = {'size' : 6}, markerscale = 100)

    
    plt.xlabel('Frequency (MHz)')
    plt.locator_params(axis="x", nbins=6)
    plt.ylabel('$|S_{21}|^2$ (dB)')
    plt.grid()
    
    #plots square magnitude fit in decibels over data
    plt.plot(freqsMHz_new, 10*np.log10(amplitudeequation(freqsMHz_new, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])),)
    
    #calculates midpt of y-axis in plot. Used for saving plot to file
    midpt = (np.min(squaremagnitudesdb) + np.max(squaremagnitudesdb)) / 2
    
    plt.savefig(pp, format='pdf')
 
   
    
    #now plot unwrapped, rotated phase data for each power
    

    plt.figure(4)
    plt.scatter(freqsMHz_new, thetafinal *  180 / np.pi , marker='.', s=0.01)
    plt.title('Phase vs Frequency')
    plt.xlabel('Frequency (MHz)')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Phase (degrees)')
    plt.plot(freqsMHz_new, phaseequation_degrees(freqsMHz_new, popt_phase[0], popt_phase[1], popt_phase[2]),linewidth = 1)
    plt.grid()
   
    plt.savefig(pp, format='pdf')

    
    #plots Qi phase with power for each power
    plt.figure(5)
    plt.scatter(power, Qi_phase, marker='.', s=100)
    plt.xlabel('Power')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Qi')
    plt.title('Qi vs Power')
    plt.grid()
    
    plt.savefig(pp, format='pdf')
    
 
    #Plots Qc phase with power for each power
    plt.figure(6)
    plt.scatter(power, Qc_phase, marker='.', s=100)
    plt.xlabel('Power')
    plt.locator_params(axis="x", nbins=9)
    plt.ylabel('Qc')
    plt.title('Qc vs Power')
    plt.grid()
    
    plt.savefig(pp, format='pdf')

    
    pp.close()
    
    
    
    
    













