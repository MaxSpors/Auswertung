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
        legend_handle.set_alpha(1)  # Setzen Sie die Opazit�t auf 1 (vollst�ndig undurchsichtig)
    ax.grid(True, color="slategray", linewidth="0.3", linestyle="--",alpha=0.4)

def feldanalyse(folder_path, Feld):
    all_subfolders = [os.path.join(folder_path, subfolder) for subfolder in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, subfolder))]
    num_subplots = len(all_subfolders)  # Anzahl der Subplots
    
    if num_subplots == 0:
        print("No subfolders found.")
        return
    
    # Extraktion der Erwartungswerte des Photopeaks
    Pos, DPos = [],[]
    E = []
    s, Ds = [],[]
    
    # Calculate the number of rows and columns for the subplots
    num_cols = int(np.ceil(np.sqrt(num_subplots)))
    num_rows = int(np.ceil(num_subplots / num_cols))
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 20), dpi=300)
    axes = axes.flatten()
    
    for ax, subfolder in zip(axes, all_subfolders):
        for root, _, files in os.walk(subfolder):
            files = [f for f in files if f.endswith('.mca')]
            files.sort(key=extract_field_strength)
            
            for i, file_path in enumerate(files):
                full_path = os.path.join(root, file_path)
                x,data, Fehler,identification = read_out(full_path)
                x = np.linspace(0, len(data), len(data))
                label = f"{identification} V/cm"
                color = plt.cm.viridis(i / len(files))
                
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
    
    # Hide any unused subplots
    for ax in axes[num_subplots:]:
        ax.axis('off')
    
    plt.suptitle(f'{Feld} Kontrollplot', fontsize=30, fontfamily='Calibri')
    plt.tight_layout()
    plt.subplots_adjust(top=0.96)
    plt.show()
    
    return Pos, DPos, E, s, Ds

def Korrektur(Array):
    filtered_Array = [x for x in Array if x is not None]

# Schritt 2: Berechnen Sie den Mittelwert der verbleibenden Werte
    mean_value = np.mean(filtered_Array)

# Schritt 3: Ersetzen Sie die None-Werte durch den berechneten Mittelwert
    Array = [x if x is not None else mean_value for x in Array]

    return Array