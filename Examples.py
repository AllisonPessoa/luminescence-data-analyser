import SpectraAnalyser as spr

#%% ### Individual spectra ###

particle1_max_power = spr.Spectrum('Example Files/od1_0_0.asc')

band1 = particle1_max_power.add_band('4F5/2', [817,824], r'$^4$F$_{5/2}$')
band2 = particle1_max_power.add_band('4F3/2', [879,886], r'$^4$F$_{3/2}$')

fig, ax = particle1_max_power.plot_spectrum()

particle1_max_power.get_area_under_bands()
#%% ### Power Dependency ###
particle_example = [
            spr.Spectrum('Example Files/od1_0_0.asc'),
            spr.Spectrum('Example Files/od1_0_2.asc'),
            spr.Spectrum('Example Files/od1_0_4.asc'),
            spr.Spectrum('Example Files/od1_0_6.asc'),
            spr.Spectrum('Example Files/od1_0_8.asc')
            ]

for spectrum in particle_example:
    spectrum.add_band('4F5/2', [817,824], r'$^4$F$_{5/2}$')
    spectrum.add_band('4F3/2', [879,886], r'$^4$F$_{3/2}$')
    
initial_power = 1e-6 #Power (W), measured by neutral density filters
filters_neutral_density = [0,0.2,0.4,0.6,0.8]
powers = [(initial_power/(10**ND)) for ND in filters_neutral_density]

power_dependence_example = spr.PowerDependence(particle_example, powers)
power_dependence_example.plot_power_law({'525 nm': [515,545], '544 nm': [545,570]},[1e-5,1e-6])

#%% ### Band Intensity Ratio ###

spectra_Amado_Batista = [
            Spectrum('Experiments/1_Amado Batista/LIR/22C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/25C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/29C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/30C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/34C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/38C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/40,8C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/43C.asc'),
            Spectrum('Experiments/1_Amado Batista/LIR/47,2C.asc'),
            ]

for i in range(len(spectra_Amado_Batista)):
    spectra_Amado_Batista[i].correct_baseline()
    #spectra_PTC01_ar[i].plot_spectrum()

temperatures = [T+273 for T in [22,25,29,30,34,38,40.8,43,47.2]]
LIR_Amado_Batista = LIR(spectra_Amado_Batista, temperatures, 'Amado Batista')

interval = {'525 nm': [517,543], '544 nm': [543,570]}
LIR_Amado_Batista.calculate(interval, stat=True) #stat = True: returns the statistics of the evaluations

fig, ax = LIR_Amado_Batista[0].plot_spectrum()
fig

