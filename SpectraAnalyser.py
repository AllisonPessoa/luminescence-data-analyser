####################################################
#           Luminescence Analysis Tools
# Python scripts used for standard analysis in Luminescence Experiments
#		Spectra analyser
# Features:
#   - Sprectral data from Monochromator + CCD camera:
#       - Power Dependency Analysis on a series of spectra for each power excitation
#       - Thermometry Analysis on a series of spectra for each temperature

# Used and tested with data from iDus 401 CCD camera
# (obtained through Andor SOLIS for Spectroscopy)

# By Allison Pessoa, Nano-Optics Laboratory.
# UFPE, Brazil.  Last update: December, 2021

# Do not share this scrip, fully or partially, without proper permission.
# allison.pessoa@ufpe.br | allisonpessoa@hotmail.com 
####################################################

import numpy as np

import scipy.integrate as integrate

import uncertainties
from uncertainties.umath import *

import matplotlib
import matplotlib.pyplot as plt

import matplotlib.font_manager

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class DataOpener():
    @staticmethod
    def from_txt(filename):
        spectrum_file = open(self.filename)
        wavelength = []
        intensity = []
      
        for line in spectrum_file.readlines():
            wavelength.append(float(line.split(',')[0]))
            intensity.append(float(line.split(',')[1]))
              
class Spectrum():
# Contains all information about each individual spectrum
    def __init__(self, file_name):
        #file_name: .asc file spectral data. File must contain two columns.
        #           first one is wavelength (decimal sep as ',')
        #           second one is intensity (counts/second, decimal sep as '.')
        #           Also accepts '.txt' files with columns separated by ','.

        self.filename = file_name
        (self.wavelength, self.intensity) = self._load_file()
        
    @classmethod
    def from_txt(cls, file_name):
        
        spectrum = cls(file_name)
        
    def _load_file(self):
        #Loads and processes the spectrum file.
        spectrum_file = open(self.filename)

        for line in spectrum_file.readlines():
            wavelength.append(float(line.split('\t')[0].replace(',','.')))
            intensity.append(float(line.split('\t')[1][:-1]))

        #Better approach if I used an Adapter (Adapter Design Pattern) - To be implemented
        elif self.filename[-4:] == '.txt':
            for line in spectrum_file.readlines():
                wavelength.append(float(line.split(',')[0]))
                intensity.append(float(line.split(',')[1]))

        spectrum_file.close()
        return (wavelength,intensity)

    def correct_baseline(self, background):
        #Subtracts the data spectrum from a given baseline spectrum.
        wvl, intensity_background = background.get_spectrum(normalized=False)
        for i in range(len(self.intensity)):
            self.intensity[i] = self.intensity[i] - intensity_background[i]

    def get_spectrum(self, normalized=False):
        #Returns a tuple with the spectral data: (wavelength, intensity)
        #If normalized is True, intensity will be normalized to values between 0 and 1
        intensity = []

        if normalized == True:
            for j in range(len(self.intensity)):
                intensity.append((self.intensity[j])/(max(self.intensity)))
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
        area = integrate.simpson(intensity_to_integrate, wavelength_to_integrate)

        return(area)

    def plot_spectrum(self, normalized=False):
        #Plots spectrum, and returns the figure and axis objects for further changes.
        #If normalized is True, intensity will be normalized to values between 0 and 1
        wavelength, intensity = self.get_spectrum(normalized = normalized)

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
# Calculates the linear dependence between excitation power and intensity
    def __init__(self, spectra_set, power_set):
        #spectra_set: list (length M) of Spectrum objects each with a different excitation power.
        #             Mininum two
        #power_set: list (length M) of the excitation power density for each spectrum
        #              (order is important).
        self.spectra_set = spectra_set
        self.power = power_set
        self.plots_cmap = plt.get_cmap("plasma")
        self.plots_CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                                     '#f781bf', '#a65628', '#984ea3',
                                     '#999999', '#e41a1c', '#dede00']
        self.plots_marker_cycle = ['o','v','*','d','^','s','>','<','g']
    
    def get_intensities(self, integration_limits_dic, normalized=False):
        #Returns a list of dictionaires including the band's name, 
        #           and the integral of the spectra in the defined limits for all powers.
        #Integration_limits_dic: Dictionary that defines two wavelength bands to compute the fittings 
        #                       Example: {'band1': [510,535], 'band2': [535,554]}
        #If normalized is True, intensity will be normalized to values between 0 and 1
        number_of_bands = len(integration_limits_dic)
        area_under_bands = []#[[p1,p2,p3,p4], [p1,p2,p3,p4], ...] ;
                           #[p1,p2,p3,p4] = same band, different powers
        
        for i in range(number_of_bands):
            area_under_bands.append([])
            for spectrum in self.spectra_set:
                area_under_bands[i].append( spectrum.get_area_under_spectrum(\
                                             list(integration_limits_dic.values())[i], \
                                                 normalized=normalized) )
        
        list_bands = []#list of dictionaries
        for i in range(len(integration_limits_dic)):
            list_bands.append(
                {'bands_name': list(integration_limits_dic.keys())[i], \
                 'intensities': area_under_bands[i]} )
            
        return list_bands
    
    def get_power_law(self, integration_limits_dic, power_limits, normalized=False):
        #Returns the linear fittings parameters A and B of the log(Intensity) vs. log(power) plot
        #                   in a defined region of power. Intensity = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
        #Integration_limits_dic: Dictionary that defines two wavelength bands to compute the fittings 
        #                       Example: {'band1': [510,535], 'band2': [535,554]}
        #Power_limits: list of lists indicating the power separation to perform the fittings
        #                       Example: [[1e2,1e3],[1e3,1e5]]
        #If normalized is True, intensity will be normalized to values between 0 and 1
        
        area_under_bands = self.get_intensities(integration_limits_dic, normalized=False)
        index_of_separations = []
        for value in power_limits:
            nearest_min_value = min(self.power, key=lambda x: abs(x-value[0]))
            nearest_max_value = min(self.power, key=lambda x: abs(x-value[1]))
            index_of_separations.append([self.power.index(nearest_min_value),\
                                         self.power.index(nearest_max_value)])
                
        list_bands = []#list of dictionaries
        for i in range(len(area_under_bands)):
            for j in range(len(power_limits)):
                log_power_to_integrate = np.log10(self.power)[index_of_separations[j][1]-1: index_of_separations[j][0]+1] #1 and 0 are inverted here because the data goes for higher to low power
                log_area_to_integrate = np.log10(area_under_bands[i]['intensities'])[index_of_separations[j][1]-1: index_of_separations[j][0]+1]
                linear_fit, cov_matrix = np.polyfit(log_power_to_integrate, log_area_to_integrate, 1,
                                cov=True) #([a[0]*x + a[1]],residual error,...)
                (A, B) = uncertainties.correlated_values([linear_fit[0], linear_fit[1]], cov_matrix)
                list_bands.append(
                        {'bands_name': list(integration_limits_dic.keys())[i], \
                        'power_region': [10**x for x in log_power_to_integrate],
                        'A': A,
                        'B': B} )
        return list_bands
        
    def plot_power_law(self, integration_limits_dic, power_limits, normalized=False):
        #Plots the integrated intensities along with the linear fittings in a log-log scale
        #Integration_limits_dic: Dictionary that defines two wavelength bands to compute the fittings 
        #                       Example: {'band1': [510,535], 'band2': [535,554]}
        #If normalized is True, intensity will be normalized to values between 0 and 1
        
        area_under_bands = self.get_intensities(integration_limits_dic, normalized=False)
        fitting_parameters = self.get_power_law(integration_limits_dic, power_limits, normalized=False)
        
        fitting_lines = []
        for i in range(len(fitting_parameters)):
            A = fitting_parameters[i]['A']
            B = fitting_parameters[i]['B']
            power_to_integrate = fitting_parameters[i]['power_region']
            fitting_lines.append([10**(x*A.n + B.n) for x in np.log10(power_to_integrate)])
                            #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)

        fig, ax = plt.subplots()
        ax.set_xlabel('Excitation Power Density (W/cm$^2$)', fontsize='x-large')
        ax.set_ylabel('Signal Integral (a.u)', fontsize='x-large')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        for i in range(len(area_under_bands)):
            intensity_to_plot = area_under_bands[i]['intensities']
            ax.plot(self.power, intensity_to_plot, color = self.plots_CB_color_cycle[i],
                                marker = self.plots_marker_cycle[i], linestyle='None', 
                                label=list(integration_limits_dic.keys())[i])
            
        for i in range(len(fitting_parameters)):
            power_to_plot = fitting_parameters[i]['power_region']
            ax.plot(power_to_plot, fitting_lines[i], color = self.plots_CB_color_cycle[i], 
                                linestyle='--', label = "Slope: {:.2f} \u00B1 {:.1}"
                                .format(fitting_parameters[i]['A'].n, fitting_parameters[i]['A'].s))
        ax.legend()
        return fig, ax

    def plot_all_spectra(self, integration_limits, normalized=False):
        #Plots all spectra in only one graph, with different collors, for comparison.
        #Integration limits is for plotting a vertical dashed line indicating the integration
        fig, ax = plt.subplots()
        ax.set_prop_cycle(color=self.plots_CB_color_cycle)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.tick_params(direction='in',which='both')
        
        for i in range(len(self.spectra_set)):
            wavelength, intensity = self.spectra_set[i].get_spectrum(normalized=normalized)
            ax.plot(wavelength,intensity, color=self.plots_cmap(len(self.spectra_set)*18-i*18),
                                            label='Power = {:.1e} W'.format(self.power[i]))
            
        for value in integration_limits:
            ax.axvline(value, ymax = 0.6, linestyle='--', color='black')
            
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
            
        return self.LIR

    def Plot_spectra_maxmin(self, integration_limits_dic, normalized = False):
        #PLots the maximum and mininum spectra, and also the integration limits used
        integration_limits = []
        for i in range(len(integration_limits_dic)):
            integration_limits += list(integration_limits_dic.values())[i]

        #Setups for spectrum plot
        fig, ax = plt.subplots(constrained_layout=True)
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        ax.set_ylabel('Intensity (counts/s)', size=8)
        ax.set_xlabel('Wavelength (nm)', size=8)
        #ax.grid(True)
        ax.tick_params(direction='in',which='both')

        wavelength, intensity_min_temp = self.spectra_set[0].get_spectrum(normalized=normalized)
        wavelength, intensity_max_temp = self.spectra_set[-1].get_spectrum(normalized=normalized)

        ax.plot(wavelength, intensity_min_temp, 'b', linewidth=1,
                            label='Temp = {:.1f} K'.format(self.temperatures[0]))
        ax.plot(wavelength, intensity_max_temp, 'r--', linewidth=1,
                            label='Temp = {:.1f} K'.format(self.temperatures[-1]))

        for value in integration_limits:
            ax.axvline(value, ymax = 0.6, linestyle='-.', linewidth = 1, color='black')
        #ax.legend()

        return fig,ax

    def plot_LIR(self):
        #Plots the dependence of the LIR with the temperature

        #Setups for decay curves plots
        fig, ax = plt.subplots(constrained_layout=True)
        ax.set_xlabel(r'1/T $(x10^{-3})$ $(\mathrm{K}^{-1})$', size=8)
        ax.set_ylabel(r'$\ln(R)$', size=8)
        ax.tick_params(direction='in',which='both')
        ax.grid()
        sec_y = ax.twinx()
        sec_y.set_ylabel(r'$S_R (\%)$', size=8, color='b')
        sec_y.tick_params(direction='in',which='both', colors='b')

        x_axis_scaled = [x*10**(3) for x in self.inverse_temperature]
        ax.scatter(x_axis_scaled, self.ln_LIR, color='000000')
        ax.plot(x_axis_scaled, self.linear_fitted_curve, color='000000')
        ax.set_title(self.particle_name, size=8)

        sens_rel = [100*x.n for x in self.relative_sens[::1]]
        sec_y.plot(x_axis_scaled, sens_rel, 'b-')

        return fig, ax
