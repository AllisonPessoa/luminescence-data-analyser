import SpectraAnalyser as spr

#%% ### Individual spectra ###

spectrum = spr.Spectrum('Example Files/od 0 + 0.0.asc')
spectrum.subtract_flat_baseline_by_region([515,535])
area = spectrum.get_area([515,535])

fig, ax = spectrum.plot_spectrum()

#%% ### Power Dependency ###
power_spectra_set = [
            spr.Spectrum('Example Files/od 0 + 0.0.asc'),
            spr.Spectrum('Example Files/od 0 + 0.2.asc'),
            spr.Spectrum('Example Files/od 0 + 0.4.asc'),
            spr.Spectrum('Example Files/od 0 + 0.6.asc'),
            spr.Spectrum('Example Files/od 0 + 0.8.asc')
            ]

filters_set = [0,0.2,0.4,0.6,0.8]
twoHeleven_power_dependence = spr.PowerDependency.load_from_ND_filters(power_spectra_set, 
                                                                       [515,535], 
                                                                       32.5e-6, 
                                                                       filters_set, 
                                                                       '$^2$H$_{11/2}$')
twoHeleven_power_dependence.plot_power_law()

#%% ### Band Intensity Ratio ###

temperature_spectra_set = [
            spr.Spectrum('Example Files/temperature 22C.asc'),
            spr.Spectrum('Example Files/temperature 30C.asc'),
            spr.Spectrum('Example Files/temperature 40C.asc'),
            spr.Spectrum('Example Files/temperature 50C.asc'),
            spr.Spectrum('Example Files/temperature 60C.asc')
            ]

temperature_set = [22,30,40,50,60]
LIR_analysis = spr.LIR_Analysis(temperature_spectra_set, [515,535], [535,545],
                        temperature_set, 'A2/A1', LIR_error=0.0022)
LIR_analysis.plot_LIR(fmt='o')
