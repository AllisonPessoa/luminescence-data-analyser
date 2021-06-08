########
#		Examples using the module
#

from SpectraAnalyser import *

########################### POWER DEPENDENCY ############################### %%
part_amado_batista = [
            Spectrum('Experiments/1_Amado Batista/PowerDependency/20_20-01-2021/OD0,0.asc'),
            Spectrum('Experiments/1_Amado Batista/PowerDependency/20_20-01-2021/OD0,1.asc'),
            Spectrum('Experiments/1_Amado Batista/PowerDependency/20_20-01-2021/OD0,2.asc'),
            Spectrum('Experiments/1_Amado Batista/PowerDependency/20_20-01-2021/OD0,4.asc'),
            Spectrum('Experiments/1_Amado Batista/PowerDependency/20_20-01-2021/OD0,6.asc')
            ]

pot_0 = 0.8*10**(-6) #Power (W), measured by neutral density filters
percent = [1, 0.74, 0.53, 0.4, 0.23]
pot = [x*pot_0 for x in percent]

filters = [0,0.1,0.2,0.4,0.6]
pot_2 = [(pot_0/(10**(i-filters[0]))) for i in filters]

power_dependence_part_amado_batista = PowerDependence(part_amado_batista, pot)
power_dependence_part_amado_batista.plot_power_dependence({'525 nm': [515,545], '544 nm': [545,570]})

########################### LIR ################################# %%

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

########################### LIFETIMES ################################# %%

decay_curves_Amado = [
            #DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 22,0C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 25,0C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 28,8C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 31,5C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 34,4C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 37,7C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 41,0C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 44,0C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 47,2C.txt'),
            DecayCurve('Experiments/1_Amado Batista/LifeTime/525nm S 50,0C.txt')
            ]

temperatures = [T+273 for T in [25,28.8,31.5,34.4,37.7,41,44,47.2,50]]
lifetime_Amado = LifeTime(decay_curves_Amado, temperatures, '525nm - Saturated')
lifetime_Amado.calculate([3450,12500],n=10)
