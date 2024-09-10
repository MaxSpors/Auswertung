# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.special import erf

def Offset(file):
    # Bereinige Readout um gemessenen Offset am PicoAmpMeter
    dataframe_currents = pd.read_csv(file, skiprows=1, sep='\t')

    # Channel 12 again
    currents = dataframe_currents.iloc[:, 28].to_numpy()
    dcurrents = dataframe_currents.iloc[:, 29].to_numpy()

    currents=currents[:-1]
    dcurrents=dcurrents[:-1]

    offset,doffset=np.mean(currents), 1/len(dcurrents)*np.sqrt(np.sum(dcurrents**2)) 
    
    return offset,doffset

def Energiekalibration(xdata,dxdata):
    # Bestimme Gerade aus den zwei Punkten der Erwartungswerte der Linien, die wir auf die Energie der Linien mappen können.
    # Bestimme Steigung als rise over run (mittlere Änderungsrate)
    m= (5.775-2.892)/(xdata[0]-xdata[1])
    dm=(5.775-2.892)/(xdata[0]-xdata[1])**2 *np.sqrt(dxdata[0]**2 +dxdata[1]**2)

    # Bestimme y-Achsenabschnitt als Offset für eine Linie
    b= 5.775-m* xdata[0]
    db=np.sqrt((xdata[0]*dm)**2 +(dxdata[0]*m)**2)
    return [m,b],[dm,db]

def calibrationfunction(x,m,b,dm,db):
    # Baue eine Gerade mit diesen Parametern und den dazugehörigen Fehlern
    return m*x+b, np.sqrt(dm**2 *x**2 + db**2)

def I_ionisation(MCAdatensatz: np.ndarray, Counts: np.ndarray, time: float, results) -> tuple:
    # Initiiere den Geradenfit zur Reskalierung der X-Achse
    Pealpos = [results.params['p1'].value, results.params['p4'].value]
    Peakerr = [results.params['p1'].stderr, results.params['p4'].stderr]
    zone_between_peaks=np.min(Counts[int(results.params['p4'].value): int(results.params['p1'].value) ])
    #zone_between_peaks=0

    Kalibparams, Kaliberr = Energiekalibration(Pealpos, Peakerr)
    EnergieMCA, dEnergieMCA = calibrationfunction(MCAdatensatz, Kalibparams[0], Kalibparams[1], Kaliberr[0], Kaliberr[1])
    for i in range(len(Counts)):
        Counts[i]=Counts[i]-zone_between_peaks*(erf((Counts[i]-results.params['p1'].value)/(results.params['p2'].value))-1)
    dCounts = np.sqrt(Counts)
    
    weighted_sum = EnergieMCA * Counts
    dweighted_sum = np.sqrt((dEnergieMCA * Counts) ** 2 + (dCounts * dEnergieMCA) ** 2)

    Summe = np.sum(weighted_sum)                                    # Egesamt
    dSumme = np.sqrt(np.sum(dweighted_sum ** 2))

    W_ion = 28.4 / 1000  # keV
    I_ioni = Summe / (W_ion* time) *1.60217662e-19
    dI_ioni = dSumme / (W_ion * time) *1.60217662e-19
    return I_ioni, dI_ioni