# -*- coding: utf-8 -*-

import numpy as np
from lmfit import Parameters, Model
import os

def Gain_best(Strom: np.ndarray, I_ionistaion: float, Offset: float, dStrom: np.ndarray, dOffset: float, dI_ionisation: float):
    # Konstante
    # Berechnung des Gains und der Unsicherheit
    Readout=np.mean(Strom)
    dReadout=np.mean(dStrom)
    Readout=Readout*10**(-12)
    dReadout=dReadout*10**(-12)
    Offset=Offset*10**(-12)
    dOffset=dOffset*10**(-12)
    Gain = np.abs((Readout - Offset) / (I_ionistaion))
    dGain = np.sqrt((dReadout / (I_ionistaion))**2 +
                    (dOffset / (I_ionistaion))**2 +
                    ((Readout - Offset) * (dI_ionisation) / (I_ionistaion)**2)**2)

    return Gain, dGain

def Gainexponential(x,A,B,C):
    return A*np.exp(B*(x-90)) +C

def get_file_names(folder_path):
    # Liste aller Dateien im angegebenen Ordner
    files = os.listdir(folder_path)
    
    # Filtere nur die Dateien (keine Verzeichnisse)
    file_names = [f for f in files if os.path.isfile(os.path.join(folder_path, f))]
    
    return file_names

def Gain_Anpassung(x: np.ndarray, y: np.ndarray, dy: np.ndarray) -> np.ndarray:
    gmodel = Model(Gainexponential)
    
    # Initialisierung der Parameter
    params = Parameters()
    params.add('A', value=100, min=80, max=600)
    params.add('B', value=0.5)
    params.add('C', value=50)

    weights=[]
    # Gewichte berechnen
    for i in range(len(dy)):
        weights.append(  1.0 / dy[i])
        
    # Fit durchführen
    result = gmodel.fit(y, params, x=x, weights=weights)
    
    return result