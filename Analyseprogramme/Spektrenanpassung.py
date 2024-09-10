# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import erf
from lmfit import Model, Parameters,fit_report


def FeSpectrum(x, p0, p1, p2, p3, p4, p5, p6, p7, p8):
    # Fitfunktion f�r die Eisenspektren
    # x: Wertebereich
    # p0,...,p8: Fitparameter f�r die einzelnen Funktionen, siehe unten

    corr1 = 6.49 / 5.89
    corr2 = 3.53 / 2.93

    # Fitfunktionen
    N = p0 * (np.exp(-1/2 * ((x - p1) / p2)**2) + 1/8.8 * np.exp(-1/2 * ((x - p1 * corr1) /( p2 * np.sqrt(corr1)))**2))               # photo peak
    N += p3 * (np.exp(-1/2 * ((x - p4) / p5)**2) + 1/8.8 * np.exp(-1/2 * ((x - p4 * corr2) / (p5 * np.sqrt(corr2)))**2))              # escape peak
    N += -p6 * (erf((x - p1)/p2) - 1)                                                                                                 # noise
    N += p7 * np.exp(-x / p8)                                                                                                         # (noise + flouresence -> Hauer)
    return N

def FWHM_calc(MCA_Array: np.ndarray, Counts_Array: np.ndarray, max_counts: float, max_index: int) -> tuple:
    #Berechnet die Halbwertsbreite und die Standardabweichung eines Peaks aus den Messdaten.
    #Wichtig: Die Messdaten m�ssen um den Offset bereits bereinigt sein.
    
    half_max = max_counts / 2

    left_points = Counts_Array[:max_index]
    right_points = Counts_Array[max_index + 1:]

    left_boundary = None
    right_boundary = None

    # Approximiere die Halbwertsbereite, indem der Wertebereich auf den Threshold gescannt wird
    for i in range(len(left_points) - 1, 0, -1):
        if left_points[i] <= half_max:
            left_boundary = max_index-(max_index-i )
            break
    for i in range(len(right_points) - 1):
        if right_points[i] <= half_max:
            right_boundary = max_index+i
            break
    if left_boundary is None or right_boundary is None:
        print("Konnte die linke oder rechte Grenze nicht finden.")
        return None

    FWHM = right_boundary - left_boundary
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))

    return FWHM, sigma

def Approximation_fit_parameters(MCA_Array: np.ndarray, Counts_array: np.ndarray) -> tuple:

    # 1. Bestimme x-Werte der Peaks �ber die feste relative Beziehung
    max_index_photo = np.argmax(Counts_array)
    max_index_escape = int(np.floor(0.481 * max_index_photo))

    # 2. Sch�tze aus dem Minimum des Zwischenbereichs die Amplitude der error-function ab
    room_between_peaks_counts = Counts_array[max_index_escape:max_index_photo]
    offset_guess = np.min(room_between_peaks_counts)

    # Bereinige die Counts um den Offset
    working_counts = Counts_array - offset_guess

    # 3. Bestimme Maximum der bereinigten Werte
    max_counts_photo = working_counts[max_index_photo]
    max_counts_escape = working_counts[max_index_escape]

    # 3. Nutze die Standardabweichung f�r die bereinigten Peaks
    FWHM_photo, sigma_photo = FWHM_calc(MCA_Array, working_counts, max_counts_photo, max_index_photo)
    FWHM_escape, sigma_escape = FWHM_calc(MCA_Array, working_counts, max_counts_escape, max_index_escape)

    # 4. Bestimme Parameter des Exponentials (setze Scan-Parameter)
    start_value =  np.max(working_counts[1:20])
    Threshold = start_value / np.sqrt(np.e)
    limit_index = max_index_escape - int( FWHM_escape)

    # F�hre den Scan durch
    time_constant = None
    for i in range(limit_index):
        if working_counts[i] <= Threshold:
            time_constant = 2 * i
            break
    if time_constant is None:
        time_constant = 10

    # 5. Approximiere Amplitude der erf()
    Amp_guess = offset_guess / 2
    time_constant=time_constant+10

    return max_counts_photo, max_index_photo, max_counts_escape, max_index_escape, sigma_photo, sigma_escape, time_constant, Amp_guess

def SpektrenAnpassung(x: np.ndarray, y: np.ndarray, dy: np.ndarray):   # Setze die Grundlagen f�r den Anpassungsprozess
    gmodel = Model(FeSpectrum)
    
    # Bestimme die Erstsch�tzer f�r die Anpassung
    max_counts_photo, max_index_photo, max_counts_escape, max_index_escape, sigma_photo, sigma_escape, time_constant, Amp_Erf = Approximation_fit_parameters(x, y)
    time_constant=time_constant+0.005

    params = Parameters()
    params.add('p0', value=max_counts_photo, min=0.6*max_counts_photo, max=1.5*max_counts_photo)
    params.add('p1', value=max_index_photo, min=max_index_photo-70, max=max_index_photo+70)
    params.add('p2', value=sigma_photo, min=0.2*sigma_photo, max=1.5*sigma_photo+90)

    params.add('p3', value=max_counts_escape, min=0.5*max_counts_escape-30, max=4*max_counts_escape+20)
    params.add('p4', value=max_index_escape, min=max_index_escape-70, max=max_index_escape+70)
    params.add('p5', value=sigma_escape, min=0.1*sigma_escape, max=2.4*sigma_escape+90)

    params.add('p6', value=Amp_Erf, min=-3*np.abs(Amp_Erf), max=3*np.abs(Amp_Erf))

    params.add('p7', value=np.max(y[0:60]), min=0.1*np.max(y[0:60]), max=1.5*np.max(y[0:60]))
    params.add('p8', value=time_constant, min=0.01*time_constant, max=4*time_constant+40)
    weights = 1.0 / dy

    result = gmodel.fit(y, params, x=x, weights=weights,method='leastsq')
    return result

