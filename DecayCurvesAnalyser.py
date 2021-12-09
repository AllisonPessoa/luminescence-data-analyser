####################################################
#               Decay Curve analyser
# Python scripts used for standard analysis in Luminescence Experiments

# Features:
#   - Decay Curves from a Time-Correlated Single Photon Counter
#       - Plot and Analyser
#       - Calculations of the Lifetime with some methods *Under Construction*
#
# Time-correlated photon counter home-assembled

# By Allison Pessoa, Nano-Optics Laboratory.
# UFPE, Brazil.
####################################################

import numpy as np

import scipy.integrate as integrate
import scipy.signal as signal
import scipy.optimize as optimize

import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class DecayCurve():
# Contains all information about each individual decay curve
    def __init__(self, file_name, sep):
        self.filename = file_name
        (self.time, self.intensity) = self._load_file(sep)

    def _load_file(self,sep):
        #Loads and processes the decay curves .txt files "time,intensity\n"
        decay_curve_file = open(self.filename)
        time = []
        intensity = []

        for line in decay_curve_file.readlines():
            time.append(float(line.split(sep)[0]))
            intensity.append(float(line.split(sep)[1]))

        decay_curve_file.close()
        return (time,intensity)
    
    def _get_sep_indexes(self, values_separation):
        #I dont know
        index_separations = []
        
        for value in values_separation:
            nearest_time_value = min(self.time, key=lambda x: abs(x-value))
            index_separations.append(self.time.index(nearest_time_value))
        
        return index_separations

    
    def apply_savgol(self, time, intensity, wl, order):
        filtered_data = signal.savgol_filter(intensity, wl, order)
        filtered_deriv = signal.savgol_filter(intensity, wl, order, deriv=1, delta=time[1]-time[0])

        return filtered_data, filtered_deriv
    
    def get_decay_curve(self, normalized=False, baseline=False, slice_data=False):
        #Returns a tuple with the decay curve data: (time, intensity)
        #If normalized is True, intensity will be normalized to values between 0 and 1
        intensity = []
        
        if baseline != False:
            indexes = self._get_sep_indexes(baseline)
            for x in self.intensity:
                intensity.append(x - np.mean(self.intensity[indexes[0]:indexes[1]]))
        else:
            intensity = self.intensity
        
        if slice_data != False:
            index_of_separations = self._get_sep_indexes(slice_data)

            intensity = intensity[index_of_separations[0]:index_of_separations[1]]
            time = self.time[index_of_separations[0]:index_of_separations[1]]
            time = [t - min(slice_data) for t in time]
        else:
            time = self.time
        
        if normalized != False:
            #mean = np.mean(intensity[normalized[0]:normalized[1]])
            intensity = [x/normalized for x in intensity]
            count_max = max(intensity)
            intensity = [(x/count_max) for x in intensity]
            
        return (time, intensity)
    
    def calculate_lifetime(self, normalized=False, baseline=False, slice_data=False, plot=True):
        
        time, intensity = self.get_decay_curve(normalized=normalized, baseline=baseline, slice_data=slice_data)
        weighted_intensity = [time[i]*intensity[i] for i in range(len(time))]
        
        weighted_area = integrate.trapezoid(weighted_intensity, time)
        norm_area = integrate.trapezoid(intensity, time)
        
        self.lifetime = weighted_area/norm_area
       
        if plot == True:
            self.plot_decay_curve(normalized=normalized, baseline=baseline, slice_data=slice_data, v_lines=[self.lifetime])
        
        return self.lifetime
    
    def plot_decay_curve(self, normalized=False, baseline=False, slice_data=False, v_lines = None, h_lines = None):
        #Plots the decay curves. Returns the matplotlib object for any alterations, if necessary.
        #If normalized is True, intensity will be normalized to values between 0 and 1

        #Setups for decay curves plots
        fig, ax = plt.subplots()
        ax.set_xlabel('Time [$\mu$s]', fontsize='x-large')
        ax.set_ylabel('Counts [a.u]', fontsize='x-large')

        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        time, intensity = self.get_decay_curve(normalized=normalized, baseline=baseline, slice_data=slice_data)
        ax.plot(time, intensity)

        if v_lines != None:
            for value in v_lines:
                ax.axvline(value, ymax = 1, linestyle='--', color='black')
        if h_lines != None:
            for value in h_lines:
                ax.axhline(value, xmax = 1, linestyle='--', color='black')

        return fig, ax

    
class LifeTime():
# Evaluates te LifeTime vs. Temperature dependence
    def __init__(self, decay_curve_set, temperature_set, particle_name):
        self.decay_curve_set = decay_curve_set
        self.temperatures = temperature_set
        self.particle_name = particle_name

    def calculate(self, time_interval, n=20, normalized=False, plot=True):
    # Calculates and plot the LT vs. T dependence
        self.lifetimes = []
        for curve in self.decay_curve_set:
            self.lifetimes.append(curve.get_lifetime(time_interval, n=n, plot=False))

        if plot == True:
            self._plot()

    def _plot(self):
        #Plots the dependence of the LIR with the temperature
        fig, ax = plt.subplots(constrained_layout=True)
        ax.scatter(self.temperatures, self.lifetimes, color='000000')

        ax.set_xlabel('Temperature (K)$', size='x-large')
        ax.set_ylabel('Lifetime ($\mu$s)', size='x-large')
        ax.tick_params(direction='in',which='both')
        ax.grid()
        ax.set_title('Lifetime dependence - Particle '+self.particle_name, size='x-large')

if __name__ == "__main__":
    print("Everything passed")