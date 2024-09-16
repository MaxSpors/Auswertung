# Funktionen zur Erstellung der Kontrollabbildungen 
import matplotlib.pyplot as plt
from Spektrenanpassung import FWHM_calc, Approximation_fit_parameters,SpektrenAnpassung,FeSpectrum
import os
import numpy as np
from Auslese import read_out, extract_field_strength

def plot_data_with_fit(ax, x, data, Fehler, label, color, alpha, fit_result=None):
    ax.errorbar(x, data, yerr=Fehler, fmt='.', label=label, color=color, alpha=alpha)
    if fit_result:
        ax.plot(x, fit_result.best_fit, '-', color=color)
    legend = ax.legend()
    for legend_handle in legend.legendHandles:
        legend_handle.set_alpha(1)  # Setzen Sie die Opazitï¿½t auf 1 (vollstï¿½ndig undurchsichtig)
    ax.grid(True, color="slategray", linewidth="0.3", linestyle="--",alpha=0.4)

import os
import numpy as np
import matplotlib.pyplot as plt



def feldanalyse(folder_path, Feld):
    # Liste aller .mca-Dateien im angegebenen Ordner
    all_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.mca')]
    
    # Sortiere die Dateien nach Größe
    all_files.sort(key=os.path.getsize)
    num_files = len(all_files)  # Anzahl der Dateien
    
    if num_files == 0:
        print("No .mca files found.")
        return
    
    # Extraktion der Erwartungswerte des Photopeaks
    Pos, DPos = [], []
    E = []
    s, Ds = [], []
    
    # Berechne die Anzahl der Zeilen und Spalten für die Subplots
    num_spectra_per_subplot = 3
    num_subplots = int(np.ceil(num_files / num_spectra_per_subplot))
    num_cols = int(np.ceil(np.sqrt(num_subplots)))
    num_rows = int(np.ceil(num_subplots / num_cols))
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 20), dpi=300)
    axes = axes.flatten()
    
    for i in range(num_subplots):
        ax = axes[i]
        start_idx = i * num_spectra_per_subplot
        end_idx = min(start_idx + num_spectra_per_subplot, num_files)
        for file_path in all_files[start_idx:end_idx]:
            x, data, Fehler, identification = read_out(file_path)
            x = np.linspace(0, len(data), len(data))
            label = f"{identification} V/cm"
            color = plt.cm.viridis(all_files.index(file_path) / len(all_files))
            
            # Fit the data
            fit_result = SpektrenAnpassung(x, data, Fehler)
            peakpos = fit_result.best_values['p1']
            peakerr = fit_result.params['p1'].stderr
            
            sigma = fit_result.best_values['p2']
            dsigma = fit_result.params['p2'].stderr

            Pos.append(peakpos)
            DPos.append(peakerr)
            E.append(identification)
            s.append(sigma)
            Ds.append(dsigma)
            
            # Plot the data and the fit
            plot_data_with_fit(ax, x, data, Fehler, label, color, 0.03, fit_result=fit_result)
    
    # Verstecke ungenutzte Subplots
    for ax in axes[num_subplots:]:
        ax.axis('off')
    
    plt.suptitle(f'{Feld} Kontrollplot', fontsize=30, fontfamily='Calibri')
    plt.tight_layout()
    plt.subplots_adjust(top=0.96)
    plt.show()
    
    return Pos, DPos, E, s, Ds


# Schritt 2: Berechnen Sie den Mittelwert der verbleibenden Werte

def Korrektur(Array):
    filtered_Array = [x for x in Array if x is not None]

# Schritt 2: Berechnen Sie den Mittelwert der verbleibenden Werte
    mean_value = np.mean(filtered_Array)

# Schritt 3: Ersetzen Sie die None-Werte durch den berechneten Mittelwert
    Array = [x if x is not None else mean_value for x in Array]

    return Array