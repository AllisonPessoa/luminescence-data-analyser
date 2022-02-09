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
import os

import scipy.integrate as integrate
import scipy

import uncertainties
from uncertainties import umath

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
            #wavelength.append(float(line.split(',')[0]))
            #intensity.append(float(line.split(',')[1]))
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
    def __init__(self, filename, normalized=False, int_time=1):
        self.filename = filename
        self.wavelength, self.intensity = Builder.get_data(filename)
        self.normalized = normalized
        self.bands = {}
    
    def loadDirectArray(wavelength, array, normalized=False, int_time=1):
        np.savetxt('file.txt', np.c_[wavelength, array], delimiter=',')
        obj = Spectrum('file.txt', normalized, int_time)
        os.remove('file.txt')
        return obj
    
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
    
    def subtract_background(self, background, integration_time=1):
        """Subtracts (and modifies) the data spectrum from a given baseline spectrum.
        background : Spectrum object
        Returns: None.
        -------"""
        wvl, intensity_background = background.get_spectrum()
        self.intensity = [(I - B)/integration_time for I, B in zip(self.intensity, intensity_background)]
    
    def subtract_flat_baseline_by_minimum(self, n_points=20):
        minimum_value = min(self.intensity)
        corresp_index = self.intensity.index(minimum_value)
        region_to_average = self.intensity[corresp_index - n_points/2 : corresp_index + n_points/2]
                       
        flat_baseline = np.mean(region_to_average)
        self.intensity = [(I - flat_baseline) for I in self.intensity]
    
    def subtract_flat_baseline_by_region(self, wvlt_interval):
        baseline_wvlt, baseline_int = self._split_data(wvlt_interval)
        flat_baseline = np.mean(baseline_int)
        self.intensity = [(I - flat_baseline) for I in self.intensity]
        
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
    
    def plot_spectrum(self, fig = None, ax = None, plot_bands = True, **kwargs):
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
        
        ax.plot(wavelength, intensity, **kwargs)

        if self.bands is not None and plot_bands is True:
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

class SpectrumGaussian():
    def __init__(self, wavelength, intensity, indexes, sigma):
        self.wavelength = wavelength
        self.intensity = intensity
        self.indexes = indexes
        self.sigma = sigma
        self.centers = {i:self.wavelength[i] for i in self.indexes}
        self.amplitudes = {i:self.intensity[i] for i in self.indexes}
    
    def gaussian(self, x, amp1,cen1,sigma1):
        return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
    
    def func(self, x, *ampl):
        intens = 0
        for i in range(len(self.indexes)):
            intens += self.gaussian(x, ampl[i], 
                                    self.centers[self.indexes[i]], 
                                    self.sigma[self.indexes[i]])
        return intens
    
    def fit_multiple_peak(self):
        popt, pcov = scipy.optimize.curve_fit(self.func, self.wavelength, 
                                              self.intensity, 
                                              p0=list(self.amplitudes.values()))
        perr = np.sqrt(np.diag(pcov))
        
        return popt, perr

    def get_all_areas(self):
        popt, perr = self.fit_multiple_peak()

        area_gaussian = {}
        for i in range(len(self.indexes)):
            fitted_y = [self.gaussian(x, popt[i], 
                                      self.centers[self.indexes[i]], 
                                      self.sigma[self.indexes[i]]) 
                        for x in self.wavelength]
            area_gaussian[self.indexes[i]] = integrate.simps(fitted_y, self.wavelength)
        
        return area_gaussian
    
    def plot_all_curves(self, fig, ax):
        popt, perr = self.fit_multiple_peak()

        for i in range(len(self.indexes)):
            fitted_y = [self.gaussian(x, popt[i], 
                                      self.centers[self.indexes[i]], 
                                      self.sigma[self.indexes[i]]) 
                        for x in self.wavelength]
            ax.plot(self.wavelength, fitted_y, '--', alpha=0.5, label=self.indexes[i])
    
        total = [self.func(x,*popt) for x in self.wavelength]
        ax.plot(self.wavelength, total, 'r')
        ax.plot(self.centers.values(), self.amplitudes.values(), 'ro')
        
    
class PowerDependency():
# Calculates the linear dependence between excitation power and intensity
    def __init__(self, band_set, power_set, label):
        #spectra_set: list (length M) of Spectrum objects each with a different excitation power.
        #             Mininum two
        #power_set: list (length M) of the excitation power density for each spectrum
        #              (order is important).
        self.label = label
        self.band_set = band_set
        self.power_set = power_set
    
    def load_from_ND_filters(band_set, initial_power, filters_set, label):
        power_set = [(initial_power/(10**ND)) for ND in filters_set]
        
        return PowerDependency(band_set, power_set, label)
        
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

    def plot_power_law(self, fig = None, ax = None, **kwargs):
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
                            label=self.label, **kwargs)

        ax.plot(self.power_set, fitting_line, 
                            linestyle='--', label = "Slope: {:.2f} \u00B1 {:.1}"
                            .format(A.n, A.s), **kwargs)
        ax.legend()
    
        return fig, ax

class LIR():
# Evaluates the LIR vs. Temperature dependence, and returns some thermometer parameters
    def __init__(self, sup_band_set, inf_band_set, temperature_set, label):
        self.sup_band_set = sup_band_set
        self.inf_band_set = inf_band_set
        self.temperature_set = temperature_set
        self.particle_name = label

    def get_LIR(self, temperature):
        assert temperature in self.temperature_set, "Input temperature value not in temperature set"
        index = self.temperature_set.index(temperature)
        sup_intensity = self.sup_band_set[index].get_area()
        inf_intensity = self.inf_band_set[index].get_area()
        LIR = (sup_intensity/inf_intensity)
        
        return LIR
    
    def get_all_LIR(self):
        LIR = []
        for T in self.temperature_set:
            LIR.append(self.get_LIR(T))
        
        return LIR
    
    def fit_Boltzmann(self):
        LIR = self.get_all_LIR()
        ln_LIR = [np.log(x) for x in LIR]
        inverse_temperature = [(1/T) for T in self.temperature_set]        
        
        linear_fit_LIR, cov = np.polyfit(inverse_temperature, ln_LIR, 1, cov=True)
    
        return linear_fit_LIR, cov
    
    def get_thermometer_params(self, print_params = False):
        
        linear_fit_LIR, cov = self.fit_Boltzmann()
        # LIR = A*exp(alpha/T); alpha = DeltaE/K_T; beta = ln(A)
        # Using uncertainties library for error hadling
        (alpha, beta) = uncertainties.correlated_values([linear_fit_LIR[0], linear_fit_LIR[1]], cov)
        Kb = 0.695034 #Boltzmann Constant, cm^-1 / K
        energy_diff = abs(alpha)*Kb

        LIR_calculated, relative_sens, deltaT = [],[],[]
        for i in range(len(self.temperature_set)):
            T = self.temperature_set[i]
            LIR_calculated.append(umath.exp(beta + alpha/T)) #'exp' function uses uncertainties.umath
            relative_sens.append(100*(abs(alpha)/(T**2)))
            deltaT.append((1/relative_sens[i].n)*(LIR_calculated[i].s/LIR_calculated[i].n)) #n -> Nominal value; s -> Uncertainty Value

        if print_params == True:
            print("alpha:   \t {:.3f}".format(alpha))
            print("beta:    \t {:.3f}".format(beta))
            print("DeltaE:  \t {:.3f}".format(energy_diff))
            print("S_r:     \t {:.3f}".format(relative_sens[0]))
            print("DeltaT:  \t {:.3f}".format(deltaT[0]))

        return alpha, beta

    def plot_LIR(self, fig = None, ax = None):
        #Plots the dependence of the LIR with the temperature
        LIR = self.get_all_LIR()
        ln_LIR = [np.log(x) for x in LIR]
        linear_fit_LIR, cov = self.fit_Boltzmann()
        inverse_temperature = [(1/T) for T in self.temperature_set]
        linear_fitted_curve = [x*linear_fit_LIR[0]+ linear_fit_LIR[1] for x in inverse_temperature]
        
        if fig is None and ax is None:
            fig, ax = plt.subplots(constrained_layout=True)
            ax.set_xlabel(r'1/T $(x10^{-3})$ $(\mathrm{K}^{-1})$', size=8)
            ax.set_ylabel(r'$\ln(R)$', size=8)
            ax.tick_params(direction='in',which='both')
            ax.grid()
            sec_y = ax.twinx()
            sec_y.set_ylabel(r'$S_R (\%)$', size=8, color='b')
            sec_y.tick_params(direction='in',which='both', colors='b')
    
        x_axis_scaled = [x*10**(3) for x in inverse_temperature]
        ax.scatter(x_axis_scaled, ln_LIR, label=self.particle_name)
        ax.plot(x_axis_scaled, linear_fitted_curve)
        
        #sens_rel = [100*x.n for x in self.relative_sens[::1]]
        #sec_y.plot(x_axis_scaled, sens_rel, 'b-')

        return fig, ax
