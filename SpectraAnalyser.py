####################################################
#           Luminescence Analysis Tools 
# Python scripts used for standard analysis in Luminescence Experiments
#		Spectrum analyser
# Features:
#   - Sprectral data from Spectrometer + CCD camera:
#       - Power Dependency Analysis on a series of spectra for each power excitation
#       - Thermometry Analysis on a series of spectra for each temperature

# Used and tested with data from iDus 401 CCD camera
# (obtained through Andor SOLIS for Spectroscopy)

# By Allison Pessoa, Nano-Optics Laboratory.
# UFPE, Brazil.  June, 2021
####################################################

import numpy as np

import scipy.integrate as integrate
import scipy.signal as signal

import uncertainties
from uncertainties.umath import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class Plots():
# Standard plot setup class

    @staticmethod
    def get_decay_baseplot():

        return fig, ax

class Spectrum():
# Contains all information about each individual spectrum
    def __init__(self, file_name):
        #file_name: .asc file spectral data. File must contain two columns.
        #           first one is wavelength (decimal sep as ',')
        #           second one is intensity (counts/second, decimal sep as '.')
        #plot: Show (True) or not the spectrum after initializing the object

        self.filename = file_name
        (self.wavelength, self.intensity) = self._load_file()

    def _load_file(self):
        #Loads and processes the spectrum file
        spectrum_file = open(self.filename)
        wavelength = []
        intensity = []

        if self.filename[-4:] == '.asc':
            for line in spectrum_file.readlines():
                wavelength.append(float(line.split('\t')[0].replace(',','.')))
                intensity.append(float(line.split('\t')[1][:-1]))

        #Better approach if I used an Adapter (Adapter Design Pattern)
        elif self.filename[-4:] == '.txt':
            for line in spectrum_file.readlines():
                wavelength.append(float(line.split(',')[0]))
                intensity.append(float(line.split(',')[1]))

        spectrum_file.close()
        return (wavelength,intensity)

    def correct_baseline(self):
        self.intensity = self.intensity - np.mean(self.intensity[0:50])

    def get_spectrum(self, normalized=False):
        #Returns a tuple with the spectral data: (wavelength, intensity)
        #If normalized is True, intensity will be normalized to values between 0 and 1
        intensity = []

        if normalized == True:
            for j in range(len(self.intensity)):
                intensity.append(self.intensity[j] + abs(min(self.intensity)))
                    #(self.intensity[j])/(max(self.intensity)-min(self.intensity)))
        else:
            intensity = self.intensity

        return (self.wavelength, intensity)

    def get_area_under_spectrum(self, integration_limits, normalized=False):
        #Returns a float with the area under the spectrum plot between integration_limits
        #Integration_limits: list of two wavelengths separating the data regions
        #                   to be integrated.
        #If normalized is True, intensity will be normalized to values between 0 and 1
        (wavelength, intensity) = self.get_spectrum(normalized=normalized)

        index_of_separations = []
        for value in integration_limits:
            nearest_wvlth_value = min(wavelength, key=lambda x: abs(x-value))
            index_of_separations.append(wavelength.index(nearest_wvlth_value))


        intensity_to_integrate = intensity[index_of_separations[0]:index_of_separations[1]]
        wavelength_to_integrate = wavelength[index_of_separations[0]:index_of_separations[1]]
        area = integrate.trapezoid(intensity_to_integrate, wavelength_to_integrate)

        return(area)

    def plot_spectrum(self, normalized=False):
        #Plots spectrum. There is a built-in plot style in class Plots. If nothing
        #               is passed as argument in base_subplot, it uses this standard
        #               You can modify the base_subplot object after, as you wish
        #If normalized is True, intensity will be normalized to values between 0 and 1
        wavelength, intensity = self.get_spectrum(normalized = normalized)

        #Setups for spectrum plot
        fig, ax = plt.subplots()
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.tick_params(direction='in',which='both')

        ax.plot(wavelength, intensity)
        ax.set_title(self.filename)

        return fig, ax

class PowerDependence():
# Calculates the linear dependency between excitation power and intensity
    def __init__(self, spectra_set, power_set):
        #spectra_set: list (length M) of Spectrum objects each with a different excitation power.
        #             Mininum two
        #power_values: list (length M) of sample excitation power in each spectrum
        #              (order is important).
        self.spectra_set = spectra_set
        self.power = power_set

    def plot_power_dependence(self, integration_limits_dic, normalized=False):
        #Evaluates the linear fitting between Power and Intensity
        #Integration_limits_dic: Dictionary that defines two wavelength bands
        #                       to perform the LIR.
        #                       Example: {'band1': [510,535], 'band2': [535,554]}
        #If normalized is True, intensity will be normalized to values between 0 and 1
        number_of_bands = len(integration_limits_dic)
        self.area_under_bands = []#[[a,b,c], [c,d,e], ...] ;
                                #[a,b,c] = same power, different intergation limits
        for i in range(number_of_bands):
            self.area_under_bands.append([])
            for spectrum in self.spectra_set:
                self.area_under_bands[i].append(spectrum.get_area_under_spectrum(list(integration_limits_dic.values())[i], normalized=normalized))

        log_areas = np.log10(self.area_under_bands)
        log_power = np.log10(self.power)

        #Setups for power-law dependency plots
        fig, ax = plt.subplots()
        ax.set_xlabel('Power Excitation [W]', fontsize='x-large')
        ax.set_ylabel('Signal Integral [a.u]', fontsize='x-large')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        marker_cycle = ['o','v','*','d']

        for i in range(number_of_bands): #for each region of interation
            linear_fit, cov_matrix = np.polyfit(log_power, log_areas[i], 1,
                            cov=True) #([a[0]*x + a[1]],residual error,...)
            fitting_line = [10**(x*linear_fit[0]+linear_fit[1]) for x in log_power]
                            #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
            uncertainties = np.sqrt(np.diag(cov_matrix))

            ax.plot(self.power, np.array(self.area_under_bands)[i], color = color_cycle[i],
                                marker = marker_cycle[i], linestyle='None', label=list(integration_limits_dic.keys())[i])
            ax.plot(self.power, fitting_line, color = color_cycle[i], linestyle='--',
                                label = "Slope: {:.1f} \u00B1 {:.1}"
                                .format(linear_fit[0], uncertainties[0])+" $[W^{-1}]$")
        ax.legend()
        return fig, ax

    def plot_all_spectra(self, integration_limits, normalized=False):
        #Plots all spectra in only one graph, with different collors, for comparison

        #Setups for spectrum plot
        fig, ax = plt.subplots()
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.tick_params(direction='in',which='both')
        for i in range(len(self.spectra_set)):
            wavelength, intensity = self.spectra_set[i].get_spectrum(normalized=normalized)
            ax.plot(wavelength,intensity, alpha=1-0.15*i,
                                            label='Power = {:.1e} W'.format(self.power[i]))

        for value in integration_limits:
            ax.axvline(value, ymax = 0.6, linestyle='--', color='black')
        #ax.legend()

        return fig, ax

class LIR():
# Evaluates the LIR vs. Temperature dependence, and returns some thermometer parameters
    def __init__(self, spectra_set, temperature_set, particle_name):
        self.spectra_set = spectra_set
        self.temperatures = temperature_set
        self.particle_name = particle_name

    def calculate(self, integration_limits_dic, normalized = False, stat = False):
        #Integration_limits_dic: Dictionary that defines two wavelength bands
        #                       to perform the LIR.
        #                       Example: {'band1': [510,535], 'band2': [535,554]}

        self.areas_band1 = []
        self.areas_band2 = []

        for spectrum in self.spectra_set:
            self.areas_band1.append(spectrum.get_area_under_spectrum(list(integration_limits_dic.values())[0], normalized=normalized))
            self.areas_band2.append(spectrum.get_area_under_spectrum(list(integration_limits_dic.values())[1], normalized=normalized))

        self.LIR = [self.areas_band1[i]/self.areas_band2[i] for i in range(len(self.temperatures))]
        self.ln_LIR = [np.log(LIR) for LIR in self.LIR]

        self.inverse_temperature = [(1/T) for T in self.temperatures]
        self.inverse_temperature_sqrd = [(1/T**2) for T in self.temperatures]

        self.linear_fit_LIR, self.cov = np.polyfit(self.inverse_temperature, self.ln_LIR, 1, cov=True)
        self.linear_fitted_curve = [x*self.linear_fit_LIR[0]+self.linear_fit_LIR[1] for x in self.inverse_temperature]

        # LIR = A*exp(alpha/T); alpha = DeltaE/K_T; beta = ln(A)
        # Using uncertainties library for error hadling
        (self.alpha, self.beta) = uncertainties.correlated_values([self.linear_fit_LIR[0], self.linear_fit_LIR[1]], self.cov)
        Kb = 0.695034 #Boltzmann Constant, cm^-1 / K
        self.energy_diff = abs(self.alpha)*Kb

        self.LIR_calculated, self.relative_sens, self.deltaT = [],[],[]
        for i in range(len(self.temperatures)):
            T = self.temperatures[i]
            self.LIR_calculated.append(exp(self.beta + self.alpha/T)) #'exp' function uses uncertainties.umath
            self.relative_sens.append((abs(self.alpha)/(T**2)))
            self.deltaT.append((1/self.relative_sens[i].n)*(self.LIR_calculated[i].s/self.LIR_calculated[i].n)) #n -> Nominal value; s -> Uncertainty Value

        if stat == True:
            print("alpha:   \t {:.3f}".format(self.alpha))
            print("beta:    \t {:.3f}".format(self.beta))
            print("DeltaE:  \t {:.3f}".format(self.energy_diff))
            print("S_r:     \t {:.3f}".format(self.relative_sens[0]))
            print("DeltaT:  \t {:.3f}".format(self.deltaT[0]))

    def Plot_spectra_maxmin(self, integration_limits_dic, normalized = False):
        #PLots the maximum and mininum spectra, and also the integration limits used
        integration_limits = []
        for i in range(len(integration_limits_dic)):
            integration_limits += list(integration_limits_dic.values())[i]

        #Setups for spectrum plot
        fig, ax = plt.subplots()
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.tick_params(direction='in',which='both')

        wavelength, intensity_min_temp = self.spectra_set[0].get_spectrum(normalized=normalized)
        wavelength, intensity_max_temp = self.spectra_set[-1].get_spectrum(normalized=normalized)

        ax.plot(wavelength, intensity_min_temp, 'b', alpha=0.5,
                            label='Temp = {:.1f} K'.format(self.temperatures[0]))
        ax.plot(wavelength, intensity_max_temp, 'r',
                            label='Temp = {:.1f} K'.format(self.temperatures[-1]))

        for value in integration_limits:
            ax.axvline(value, ymax = 0.6, linestyle='--', color='black')
        ax.legend()

        return fig,ax

    def plot_LIR(self):
        #Plots the dependence of the LIR with the temperature

        #Setups for decay curves plots
        fig, ax = plt.subplots(constrained_layout=True)
        ax.set_xlabel(r'1/T $(x10^{-3})$ $(\mathrm{K}^{-1})$', size='x-large')
        ax.set_ylabel(r'$\ln(R)$', size='x-large')
        ax.tick_params(direction='in',which='both')
        ax.grid()
        sec_y = ax.twinx()
        sec_y.set_ylabel(r'$S_R (\%)$', size='x-large', color='b')
        sec_y.tick_params(direction='in',which='both', colors='b')

        x_axis_scaled = [x*10**(3) for x in self.inverse_temperature]
        ax.scatter(x_axis_scaled, self.ln_LIR, color='000000')
        ax.plot(x_axis_scaled, self.linear_fitted_curve, color='000000')
        ax.set_title('Nanothermometer - '+self.particle_name, size='x-large')

        sens_rel = [100*x.n for x in self.relative_sens[::1]]
        sec_y.plot(x_axis_scaled, sens_rel, 'b-')

        return fig, ax
