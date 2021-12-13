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
#from uncertainties.umath import *

import matplotlib
import matplotlib.pyplot as plt

import matplotlib.font_manager

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class DataLoader():
    @staticmethod
    def from_txt(file):
        wavelength, intensity = [], []
        for line in file.readlines():
            wavelength.append(float(line.split(',')[0]))
            intensity.append(float(line.split(',')[1]))
            
        return (wavelength, intensity)
    
    @staticmethod
    def from_asc(file):
        wavelength, intensity = [], []

        for line in file.readlines():
            wavelength.append(float(line.split('\t')[0].replace(',','.')))
            intensity.append(float(line.split('\t')[1][:-1]))
            
        return (wavelength, intensity)
    
    @staticmethod
    def from_fits(file):
        wavelength, intensity = [], []
    
        for line in file.readlines():
            wavelength.append(float(line.split('\t')[0].replace(',','.')))
            intensity.append(float(line.split('\t')[1][:-1]))
            
        return (wavelength, intensity)

class Builder():
    _available_loaders = {'txt': DataLoader.from_txt,
                          'asc': DataLoader.from_asc,
                          'fits': DataLoader.from_fits}
    
    def __init__(self, filename):
        self.file = open(filename)
        self.get_format(filename)
        
    def get_format(self, filename):
        self.file_format =  filename[filename.rfind('.') + 1 :]
    
    def extract_data(self):
        x, y = self._available_loaders[self.file_format](self.file)
        return x, y
    
    def close(self):
        self.file.close()
        
    @classmethod
    def get_data(cls, filename):
        builder = cls(filename)
        x, y = builder.extract_data()
        builder.close()
        
        return x, y

class Spectrum():
    def __init__(self, filename, normalized=False, aq_time=""):
        self.filename = filename
        self.wavelength, self.intensity = Builder.get_data(filename)
        self.normalized = normalized
        self.bands = {}

    def _split_data(self, wavelength_values):
        """Get the index of the nearest wavelength from the argument in the list.
        - wavelength_values : list of floats
        Returns: List of indexes of the same size as the argument
        """
        assert len(wavelength_values) == 2
        indexes = []
        (wavelength, intensity) = self.get_spectrum()
        
        for value in wavelength_values:
            nearest_wvlth_value = min(wavelength, key=lambda x: abs(x-value))
            indexes.append(wavelength.index(nearest_wvlth_value))
        
        band_wavelength = wavelength[indexes[0] : indexes[1]]
        band_intensity = intensity[indexes[0] : indexes[1]]
        
        return band_wavelength, band_intensity
    
    def correct_baseline(self, background):
        """Subtracts (and modifies) the data spectrum from a given baseline spectrum.
        background : Spectrum object
        Returns: None.
        -------"""
        wvl, intensity_background = background.get_spectrum()
        self.intensity = [(I - B) for I, B in zip(self.intensity, intensity_background)]
    
    def get_spectrum(self):
        """Returns a tuple with the spectral data: (wavelength, intensity)
        - normalized : Bool, optional
                The default is False. If normalized is True, 
                intensity will be normalized to values between 0 and 1
        Returns 
        Wavelengths : List
        Intensity : List
        """
        if self.normalized == True:
            intensity = [value/max(self.intensity) for value in self.intensity]
        else:
            intensity = self.intensity
            
        return (self.wavelength, intensity)
    
    def add_band(self, name, limits, label=''):
        
        x_data, y_data = self._split_data(limits)
        band = SpectrumBand(x_data, y_data, name, label)
        self.bands[name] = band
        return band
        
    def get_area_under_bands(self):
        area = {}
        for keys in self.bands:
            area[keys] = self.bands[keys].get_area()
        
        return area
    
    def plot_spectrum(self, fig = None, ax = None):
        """Plots spectrum, and returns the figure and axis objects for further changes.
        - normalized : Bool, optional. If True, intensity will be normalized 
                to values between 0 and 1
            DESCRIPTION. The default is False.
        """
        wavelength, intensity = self.get_spectrum()
        
        if fig is None and ax is None:
            fig, ax = plt.subplots()
            ax.set_ylabel('Intensity (a.u)', size='x-large')
            ax.set_xlabel('Wavelength (nm)', size='x-large')
            ax.tick_params(direction='in',which='both')
            ax.set_title(self.filename)
        
        ax.plot(wavelength, intensity, color='k')

        if self.bands is not None:
            for keys in self.bands:
                self.bands[keys].add_plot_band(fig, ax)
        
        return fig, ax

class SpectrumBand(Spectrum):
    def __init__(self, wavelength, intensity, name, label):
        self.wavelength, self.intensity = wavelength, intensity
        self.name = name
        self.label = label
    
    def add_plot_band(self, fig, ax):
        ax.axvline(self.wavelength[0], color='k', linestyle='--')
        ax.axvline(self.wavelength[-1], color='k', linestyle='--')
        ax.fill_between(self.wavelength, self.intensity, label=self.label)
        ax.legend()
    
    def get_area(self):
        """Returns a float with the area under the spectrum plot between integration_limits
        - integration_limits : list of two floats
                Two wavelengths separating the data regions to integration
        """
        area = integrate.simpson(self.intensity, self.wavelength)
        return area

class PowerDependence():
# Calculates the linear dependence between excitation power and intensity
    def __init__(self, band_set, power_set, label):
        #spectra_set: list (length M) of Spectrum objects each with a different excitation power.
        #             Mininum two
        #power_set: list (length M) of the excitation power density for each spectrum
        #              (order is important).
        self.label = label
        self.band_set = band_set
        self.power_set = power_set

    def get_area_under_bands(self):
        area_under_bands = [band.get_area() for band in self.band_set]
        return area_under_bands
    
    def get_power_law_coeffs(self):
        #Returns the linear fittings parameters A and B of the log(Intensity) vs. log(power) plot
        #                   in a defined region of power. Intensity = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
        
        area_under_bands = self.get_area_under_bands()
        log_power = np.log10(self.power_set)
        log_area = np.log10(area_under_bands)
        linear_fit, cov_matrix = np.polyfit(log_power, log_area, 1,
                        cov=True) #([a[0]*x + a[1]],residual error,...)
        (A, B) = uncertainties.correlated_values([linear_fit[0], linear_fit[1]], cov_matrix)
        
        return A, B
    
    def get_fitting_line(self):
        A, B = self.get_power_law_coeffs()
        fitting_line = [10**(x*A.n + B.n) for x in np.log10(self.power_set)]
                            #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
        return fitting_line

    def plot_power_law(self, fig = None, ax = None):
        #Plots the integrated intensities along with the linear fittings in a log-log scale
        area_under_bands = self.get_area_under_bands()
        fitting_line  = self.get_fitting_line()
        A, B = self.get_power_law_coeffs()
        
        if fig is None and ax is None:
            fig, ax = plt.subplots()        
            ax.set_xlabel('Excitation Power Density (W/cm$^2$)', fontsize='x-large')
            ax.set_ylabel('Signal Integral (a.u)', fontsize='x-large')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(which='both')    
            ax.tick_params(direction='in',which='both')


        ax.plot(self.power_set, area_under_bands,
                            marker = 'o', linestyle='None', 
                            label=self.label)

        ax.plot(self.power_set, fitting_line, 
                            linestyle='--', label = "Slope: {:.2f} \u00B1 {:.1}"
                            .format(A.n, A.s))
        ax.legend()
    
        return fig, ax

    # def plot_all_spectra(self, integration_limits, normalized=False):
    #     #Plots all spectra in only one graph, with different collors, for comparison.
    #     #Integration limits is for plotting a vertical dashed line indicating the integration
    #     fig, ax = plt.subplots()
    #     ax.set_prop_cycle(color=self.plots_CB_color_cycle)
    #     ax.set_ylabel('Intensity (a.u)', size='x-large')
    #     ax.set_xlabel('Wavelength (nm)', size='x-large')
    #     ax.grid(True)
    #     ax.tick_params(direction='in',which='both')
        
    #     for i in range(len(self.spectra_set)):
    #         wavelength, intensity = self.spectra_set[i].get_spectrum(normalized=normalized)
    #         ax.plot(wavelength,intensity, color=self.plots_cmap(len(self.spectra_set)*18-i*18),
    #                                         label='Power = {:.1e} W'.format(self.power[i]))
            
    #     for value in integration_limits:
    #         ax.axvline(value, ymax = 0.6, linestyle='--', color='black')
            
    #     return fig, ax

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