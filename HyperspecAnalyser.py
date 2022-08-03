import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append(r'/home/allison/Dropbox/4_Projetos/Projetos em andamento/Python Modules/Luminescence Data Analyser')
import SpectraAnalyser as spr

class Hyperspec():
    def __init__(self, fileName):
        data = np.load(fileName, allow_pickle=True)
        wavelength = data[0]
        hyperData = data[1]

        self.dim = len(hyperData)
        self.hyperSpectrum = []
        
        for i in range(self.dim):
            self.hyperSpectrum.append([])
            for j in range(self.dim):
                self.hyperSpectrum[i].append(spr.Spectrum.loadDirectArray(wavelength, hyperData[i,j]))

    def removeBackground(self, bgd_region):
        for i in range(self.dim):
            for j in range(self.dim):
                self.hyperSpectrum[i][j].subtract_flat_baseline_by_region(bgd_region)
    
    def transpose(self):
        self.hyperSpectrum = np.transpose(self.hyperSpectrum)
        
    def correctHisteresis(self, hist_pixels):
        mask = []
        
        for i in range(self.dim):
            mask.append([])
            for j in range(self.dim):
                if (i%2 == 0) and (j < self.dim - hist_pixels):
                    mask[i].append(self.hyperSpectrum[i][j+hist_pixels])
                else:
                    mask[i].append(self.hyperSpectrum[i][j])
    
        self.hyperSpectrum = np.copy(mask)

    def center_by_maximum_intensity(self, band_region):

        mask = self.getIntensityMap('overall', band_region)
        maximum  =  np.amax(mask)
        loc = np.where(mask == maximum)
        
        mask2 = []
        
        for i in range(self.dim):
            mask2.append([])
            for j in range(self.dim):
                try:
                    mask2[i].append( self.hyperSpectrum[i - (int(self.dim/2)-int(loc[0][0]))][j - (int(self.dim/2)-int(loc[1][0]))] )
                except:
                    mask2[i].append(self.hyperSpectrum[-1][-1])
         
        self.hyperSpectrum = np.copy(mask2)
    
    def getIntensityMap(self, band_name, band_region):
        mask = []

        for i in range(self.dim):
            mask.append([])
            for j in range(self.dim):
                area = self.hyperSpectrum[i][j].add_band(band_name, band_region).get_area()
                mask[i].append(area)
    
        return mask

    def getLIRMap(self, inf_band_name, inf_band_region, sup_band_name, sup_band_region, min_criteria):
        LIR = []

        for i in range(self.dim):
            LIR.append([])
            for j in range(self.dim):
                area_sup = self.hyperSpectrum[i][j].add_band(sup_band_name, sup_band_region).get_area()
                area_inf = self.hyperSpectrum[i][j].add_band(inf_band_name, inf_band_region).get_area()

                if (area_sup + area_inf) > min_criteria:
                    LIR[i].append(area_sup/area_inf)
                else:
                    LIR[i].append(0)

        return LIR


    def getTempMap(self, inf_band_name, inf_band_region, sup_band_name, sup_band_region, alpha, beta):
        Temp = []

        for i in range(self.dim):
            Temp.append([])
            for j in range(self.dim):
                area_sup = self.hyperSpectrum[i][j].add_band(sup_band_name, sup_band_region).get_area()
                area_inf = self.hyperSpectrum[i][j].add_band(inf_band_name, inf_band_region).get_area()

                R = (area_sup/area_inf)
                T = abs(alpha)*1/(beta - np.log(R))
                Temp[i].append(T)
            
        return Temp