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
from scipy.optimize import curve_fit
import scipy

import uncertainties
from uncertainties import umath

import matplotlib
import matplotlib.pyplot as plt

from prettytable import PrettyTable

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class DataLoader():
    """Implements the routines to separate the spectral data for each given data format"""
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
            try:
                wavelength.append(float(line.split(',')[0]))
                intensity.append(float(line.split(',')[1]))
            except:
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
    """Builder to choose the correct DataLoader function"""
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
    def __init__(self, filename):
        """Loads the spectrum through the file containing the data.
        The class DataLoader implements the specific routine for each desired data format
        filename: (string) e.g. file.txt or file.asc
        """
        self.filename = filename
        if filename:
            self.wavelength, self.intensity = Builder.get_data(filename)
        self.bands = {}

    def loadDirectArray(wavelength, intensity):
        """Creates a Spectrum instance direct from lists
        wavelength, intensity : lists of the same size.
        Returns: The Spectrum object
        """
        obj = Spectrum('')
        obj.wavelength = wavelength
        obj.intensity = intensity

        return obj

    def _split_data(self, wavelength_values):
        """Get the index of the nearest wavelength from the argument in the list.
        - wavelength_values : list of two floats (wavelength)
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
        """Subtracts the spectral data from a given baseline spectrum.
        background : Spectrum object
        Returns: None.
        -------"""
        wvl, intensity_background = background.get_spectrum()
        self.intensity = [(I - B)/integration_time for I, B in zip(self.intensity, intensity_background)]

    def subtract_flat_baseline_by_minimum(self, n_points=20):
        """Subtracts the spectral data from a constant average value n_points around the minimum intensity value.
        n_points : integer
        Returns: (Float) the value subtracted
        -------"""
        minimum_value = min(self.intensity)
        corresp_index = self.intensity.index(minimum_value)
        region_to_average = self.intensity[corresp_index - n_points/2 : corresp_index + n_points/2]

        flat_baseline = np.mean(region_to_average)
        self.intensity = [(I - flat_baseline) for I in self.intensity]
        return flat_baseline

    def subtract_flat_baseline_by_region(self, wvlt_interval):
        """Subtracts the spectral data from the average intensity value in wvlt_interval
        wvlt_interval : list of two wavelengths values
        Returns: (Float) the value subtracted
        -------"""
        baseline_wvlt, baseline_int = self._split_data(wvlt_interval)
        flat_baseline = np.mean(baseline_int)
        self.intensity = [(I - flat_baseline) for I in self.intensity]
        return flat_baseline

    def normalize_spectrum(self, wvlt_interval):
        """Divide the intensities by the maximum value in wvlt_interval
        wvlt_interval : list of two wavelengths values
        Returns: (Float) the maximum value
        -------"""
        max_wvlt, max_int = self._split_data(wvlt_interval)
        maximum_counts = np.max(max_int)
        norm_intensity = [(I/maximum_counts) for I in self.intensity]
        return Spectrum.loadDirectArray(self.wavelength, norm_intensity)

    def get_spectrum(self):
        """Returns the spectral data: (wavelenght list, intensity list)
        Returns: (list, list)
        -------"""
        return (self.wavelength, self.intensity)

    def get_area(self, wvlt_interval):
        """Returns: (float) the area under the defined limits
        limits : list of two wavelengths values
        -------"""
        x_data, y_data = self._split_data(wvlt_interval)
        area = integrate.simpson(y=y_data, x=x_data)
        return area

    def get_spectrum_in_energy_space(self):
        """Returns: (float) the wavelength barycenter under the defined limits
        limits : list of two wavelengths values
        -------"""
        hc = 1240 #eV*nm
        energy_interval = [hc/wavelength for wavelength in self.wavelength]
        jacobian = [hc/energy**2 for energy in energy_interval]
        energy_spectrum = np.multiply(self.intensity, jacobian)
        
        return Spectrum.loadDirectArray(energy_interval, energy_spectrum)
            
    def get_barycenter(self, wvlt_interval):
        """Returns: (float) the barycenter under the defined limits
        limits : list of two wavelengths values
        -------"""
        x_data, y_data = self._split_data(wvlt_interval)
        weighted_area = integrate.simpson(y=np.multiply(y_data, x_data), x=x_data)
        norm_area = integrate.simpson(y=y_data, x=x_data)
        barycenter = weighted_area/norm_area
        return barycenter
        
    def plot_spectrum(self, fig = None, ax = None, **kwargs):
        """Plots the spectral data intensity vs. wavelength. If fig and ax are passed as arguments, simply adds this datapoints to the plot. If not, creates a Fig and Axis objects
        fig, ax : matplotlib objects Figure and Axis
        Returns: (fig, ax) created or the same passed as argument"""
        wavelength, intensity = self.get_spectrum()

        if fig is None and ax is None:
            fig, ax = plt.subplots()
            ax.set_ylabel('Intensity (a.u)', size='x-large')
            ax.set_xlabel('Wavelength (nm)', size='x-large')
            ax.tick_params(direction='in',which='both')
            ax.set_title(self.filename)

        ax.plot(wavelength, intensity, **kwargs)

        return fig, ax

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


class FittingMethods():
    """Implements the methods to fit the power-law data"""
    @staticmethod
    def func_power_law(x, A, n):
        """'simple' power law: I = A*P^n, P is the excitation power"""
        return A*x**(n)

    @staticmethod
    def func_power_saturation(x, a, w):
        """second-order-process saturation power law: I = A*(P^2/(1+w*P)), P is the excitation power"""
        return ( a*(x**2/(1 + w*x)) )

class PowerDependency():
    """Calculates the dependence between excitation power and emission intensity"""
    def __init__(self, spectra_set, power_set, wvlt_interval, label, method = FittingMethods.func_power_law):
        """Loads the data.
        spectra_set: list with the Spectrum objects
        power_set: list with the corresponding powers (float)
        wvlt_interval : list of two wavelengths values (float)
        label: String for legend
        method: power-law function by default. It can be passed any function in FittingMethods class
        """
        self.label = label
        self.spectra_set = spectra_set
        self.power_set = power_set
        self.wvlt_interval = wvlt_interval
        self.method = method

    def load_from_ND_filters(spectra_set, wvlt_interval, initial_power, filters_set, label, method = FittingMethods.func_power_law):
        """Creates a PowerDependency instance direct from the list of Neutral Density (Optical Density - OD) filters
        spectra_set: list with the Spectrum objects
        wvlt_interval : list of two wavelengths values (float)
        initial_power: the power corresponding to the first spectrum (float)
        filters_set : list with the corresponding OD filters (float)
        Returns: The PowerDependency object
        """
        power_set = [(initial_power/(10**ND)) for ND in filters_set]
        return PowerDependency(spectra_set, power_set, wvlt_interval, label, method)

    def get_area_under_bands(self):
        """Get the areas under all spectra
        Returns: a list with the corresponding area for each spectrum
        """
        area_under_bands = [spectrum.get_area(self.wvlt_interval) for spectrum in self.spectra_set]
        return area_under_bands

    def get_fitting_params(self, p0=[1e7,2]):
        """ Returns the fitting parameters, depending on the method used.
        Returns: tuple with the params values with their uncertainties
        """
        area_under_bands = self.get_area_under_bands()
        popt, pcov = curve_fit(self.method, self.power_set, area_under_bands, p0, maxfev=5000)
        params = uncertainties.correlated_values(popt, pcov)
        return params

    def get_fitting_line(self):
        """ Returns the fitting lines (list of intensities), depending on the method used.
        """
        params = [param.n for param in self.get_fitting_params()]
        line = [self.method(power, *params) for power in self.power_set]
        return line

    def plot_power_law(self, fig = None, ax = None, **kwargs):
        """Plots the integrated intensities vs. power in a log-log scale. If fig and ax are passed as arguments, simply adds this datapoints to the plot. If not, creates a Fig and Axis objects
        fig, ax : matplotlib objects Figure and Axis
        """
        area_under_bands = self.get_area_under_bands()
        fitting_line  = self.get_fitting_line()

        if fig is None and ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel('Excitation Power Density (W/cm$^2$)', fontsize='x-large')
            ax.set_ylabel('Signal Integral (a.u)', fontsize='x-large')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(which='both')
            ax.tick_params(direction='in',which='both')

        ax.scatter(self.power_set, area_under_bands, s=5, label=self.label, **kwargs)
        ax.plot(self.power_set, fitting_line, linestyle='--', markersize=5, **kwargs)
        ax.legend()

        return fig, ax

class LIR_Analysis():
    """Evaluates the LIR vs. Temperature dependence, and returns some thermometer parameters
    """
    def __init__(self, spectra_set, sup_wvlt_interval, inf_wvlt_interval, temperature_set, label, LIR_error = 0.22):
        self.spectra_set = spectra_set
        self.sup_wvlt_interval = sup_wvlt_interval
        self.inf_wvlt_interval = inf_wvlt_interval
        self.temperature_set = temperature_set
        self.particle_name = label
        self.LIR_error = LIR_error

    def get_LIR(self, temperature):
        """Returns the LIR value for the given temperature, which has to be in the temperature set
        """
        assert temperature in self.temperature_set, "Input temperature value not in temperature set"
        index = self.temperature_set.index(temperature)
        sup_intensity = self.spectra_set[index].get_area(self.sup_wvlt_interval)
        inf_intensity = self.spectra_set[index].get_area(self.inf_wvlt_interval)
        LIR = uncertainties.ufloat(sup_intensity/inf_intensity, self.LIR_error) #nominal value and standad deviation

        return LIR

    def get_all_LIR(self):
        """Returns the LIR set for all spectra in spectra_set
        """
        LIR = [self.get_LIR(T) for T in self.temperature_set]
        return LIR

    @staticmethod
    def fit_linear(x_data, y_data, sigma):
        """Direct implementation of the linear ftting algorithm in which each y_data has its own uncertainty. y_data = a_0 + a_1*x_data
        x_data : list of floats
        y_data : list of floats
        sigma : list with uncertainties in y_data for each point
        Returns the fitting parameters and the covariance matrix ([float, float], cov)
        """
        inv_sigma_sqr = [1/(s)**2 for s in sigma]
        u_00 = np.sum(inv_sigma_sqr)
        u_11 = np.sum(np.multiply(np.power(x_data,2), inv_sigma_sqr))
        u_10 = np.sum(np.multiply(x_data, inv_sigma_sqr))
        v_0 = np.sum(np.multiply(y_data, inv_sigma_sqr))
        v_1 = np.sum(np.multiply(np.multiply(y_data, x_data), inv_sigma_sqr))

        Delta = u_00*u_11 - np.power(u_10,2)
        a_0 = (u_11*v_0 - u_10*v_1)/Delta
        a_1 = (-u_10*v_0 + u_00*v_1)/Delta
        cov = np.multiply([[u_11, (-1)*u_10],[(-1)*u_10, u_00]], 1/Delta)

        return [a_0, a_1], cov

    def fit_Boltzmann(self):
        """Adjust the data to be fitted by linear relation
        LIR = A*exp(alpha/T); alpha = DeltaE/K_T; beta = ln(A)
        Returns: ([float, float], matrix) parameters and covariance matrix
        """
        LIR = self.get_all_LIR()
        ln_LIR = [umath.log(x) for x in LIR]
        inverse_temperature = [(1/T) for T in self.temperature_set]

        popt, pcov = LIR_Analysis.fit_linear(inverse_temperature, [x.n for x in ln_LIR], [x.s for x in ln_LIR])
        return popt, pcov

    def get_thermometer_params(self, print_params = False):
        """Returns the thermometer parameters alpha and beta. If print_params is True, print a table with the complete thermometer information
        Returns: (float, float) alpha and beta parameters
        """
        popt, pcov = self.fit_Boltzmann()
        (beta, alpha) = uncertainties.correlated_values([popt[0], popt[1]], pcov) # Using uncertainties library for error hadling
        Kb = 0.695034 #Boltzmann Constant, cm^-1 / K
        energy_diff = abs(alpha)*Kb

        LIR = self.get_all_LIR()
        T_calculated = [(alpha/(umath.log(R)-beta)) for R in LIR]
        relative_sens = [(100*(abs(alpha)/(T**2))) for T in self.temperature_set]

        #LIR_calculated.append(umath.exp(beta + alpha/T)) #'exp' function uses uncertainties.umath
        #deltaT.append((100/relative_sens[i].n)*(LIR_calculated[i].s/LIR_calculated[i].n)) #n -> Nominal value; s -> Uncertainty Value

        if print_params == True:
            x = PrettyTable()
            x.field_names = ["PTC", "Alpha", "Beta", "Energy Diff (cm^-1)","Rel. Sens (%)","Init Temp. Calc. (K)"]
            x.add_row([self.particle_name, alpha, beta, energy_diff, relative_sens[0], T_calculated[0]])
            print(x)

        return alpha, beta

    def plot_LIR(self, fig = None, ax = None, fmt=None, direct=False, **kwargs):
        """Plots the dependence of the LIR with the temperature in a ln(LIR) vs. 1/T plot. If fig and ax are passed as arguments, simply adds this datapoints to the plot. If not, creates a Fig and Axis objects
        fig, ax : matplotlib objects Figure and Axis
        Returns: (fig, ax) created or the same passed as argument
        """
        LIR = self.get_all_LIR()
        ln_LIR = [umath.log(x) for x in LIR]
        popt, pcov = self.fit_Boltzmann()
        inverse_temperature = [(1/T) for T in self.temperature_set]
        linear_fitted_curve = [x*popt[1]+ popt[0] for x in inverse_temperature]

        if direct is False:
            if fig is None and ax is None:
                fig, ax = plt.subplots(constrained_layout=True)
                ax.set_xlabel(r'1/T $(x10^{-3})$ $(\mathrm{K}^{-1})$', size=8)
                ax.set_ylabel(r'$\ln(R)$', size=8)
                ax.tick_params(direction='in',which='both')
                ax.grid()

            x_axis_scaled = [x*10**(3) for x in inverse_temperature]
            ax.errorbar(x=x_axis_scaled, y=[x.n for x in ln_LIR], yerr=[x.s for x in ln_LIR], label=self.particle_name, **kwargs)
            ax.plot(x_axis_scaled, linear_fitted_curve, **kwargs)
        else:
            if fig is None and ax is None:
                fig, ax = plt.subplots(constrained_layout=True)
                ax.set_xlabel(r'Temperature (K)', size=8)
                ax.set_ylabel(r'LIR', size=8)
                ax.tick_params(direction='in',which='both')
                ax.grid()
            x_axis_scaled = [x*10**(3) for x in inverse_temperature]
            ax.errorbar(x=self.temperature_set, y=[x.n for x in LIR], yerr=[x.s for x in LIR], label=self.particle_name, **kwargs)

        return fig, ax
