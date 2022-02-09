import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append(r'/home/allison/Dropbox/4_Projetos/Projetos em andamento/Python Modules/Luminescence Data Analyser')
import SpectraAnalyser as spr

class Hyperspec():
    def __init__(self, fileName):
        self._load(fileName)
    
    def _load(self, fileName):
        data = np.load(fileName, allow_pickle=True)
        wavelength = data[0]
        hyperData = data[1]

        self.dim = len(hyperData)
        self.hyperSpectrum = []
        
        for i in range(self.dim):
            self.hyperSpectrum.append([])
            for j in range(self.dim):
                self.hyperSpectrum[i].append(spr.Spectrum.loadDirectArray(wavelength, hyperData[i,j]))
                self.hyperSpectrum[i][j].subtract_flat_baseline_by_region([588,590])

        del data
    
    def getIntensityMap(self, band_name, band_region):
        mask = []

        for i in range(self.dim):
            mask.append([])
            for j in range(self.dim):
                area = self.hyperSpectrum[i][j].add_band(band_name, band_region).get_area()
                mask[i].append(area)
        mask = np.transpose(mask)
        
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
        LIR = np.transpose(LIR)
        
        return LIR
    
    def correctHisteresis(self, mask, pixels):
        mask2 = np.zeros((self.dim,self.dim))
        self.hist_pixels = pixels
        for i in range(self.dim):
            for j in range(self.dim-pixels):
                if i%2 == 0:
                    mask2[i][j] = mask[i][j+pixels]
                else:
                    mask2[i][j] = mask[i][j]
                    
        return mask2
    
    def getSpectrum(self, k, l):
        spec_mask = np.transpose(self.hyperSpectrum)
        spec_mask2 = []
    	
        for i in range(self.dim):
            spec_mask2.append([])
            for j in range(self.dim-self.hist_pixels):
                if i%2 == 0:
                    spec_mask2[i].append(spec_mask[i][j+self.hist_pixels])
                else:
                    spec_mask2[i].append(spec_mask[i][j])
                    
        return spec_mask2[k][l]