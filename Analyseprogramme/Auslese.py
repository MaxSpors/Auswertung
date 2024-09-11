# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import re
import os
from datetime import datetime
# Wir unterscheiden 2 F�lle .txt files und .mca files

def extract_field_strength(filename):
    # Extrahiert die Feldstaerke aus dem Dateinamen um die Plots zu benennen und zu ordnen
    patterns = [r'Drift_(\d+)Vcm', r'TF2_(\d+)Vcm', r'TF1_(\d+)Vcm',r'GEM1_(\d+)V', r'GEM2_(\d+)V',r'Ind_(\d+) Vcm']

    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            return int(match.group(1))

    return 0

def extract_timestamp(filename):
    # Extrahiert den Timestamp aus dem Dateinamen
    timestamp_pattern = r'(\d{8}_\d{6})'
    match = re.search(timestamp_pattern, filename)
    if match:
        timestamp_str = match.group(1)
        # Konvertiert den extrahierten String in ein datetime-Objekt
        timestamp = datetime.strptime(timestamp_str, '%Y%m%d_%H%M%S')
        return timestamp
    return None


def read_out(file_path):
    file_name = os.path.basename(file_path)
    if file_path.endswith('.mca'):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
            # Auslesebereich einstellen
            data_start_index = lines.index('<<DATA>>\n') + 1
            data_end_index = lines.index('<<END>>\n')
            data_lines = lines[data_start_index:data_end_index]
            data = [int(line.strip()) for line in data_lines]
            
            df = pd.DataFrame(data, columns=['Counts'])
            counts = df.to_numpy().flatten()
            dCounts=np.sqrt(counts)
            identification=extract_field_strength(file_name)
    else:
        dataframe_spectrum = pd.read_csv(file_path)
        counts = dataframe_spectrum.iloc[:, 0].to_numpy()
        counts=counts[2:]
        dCounts = np.sqrt(counts)
        identification=extract_timestamp(file_name)

    mask = counts >= 3
    #MCA = MCA[mask]
    dCounts = dCounts[mask]
    counts = counts[mask]
    MCA = np.arange(0, len(counts))

    return MCA, counts, dCounts, identification


def process_photopeak_and_field_strength(Photopeakposition, Feldstaerke, target_value=2400):
    # Umwandeln der Listen in NumPy-Arrays
    Photopeakposition = np.array(Photopeakposition)
    Feldstaerke = np.array(Feldstaerke)

    # Finde die Indizes, bei denen Feldstaerke den target_value hat
    indices = np.where(Feldstaerke == target_value)

    # Extrahiere die entsprechenden Werte aus Photopeakposition
    values_to_average = Photopeakposition[indices]

    # Berechne den Mittelwert
    average_value = np.mean(values_to_average)

    # Ersetze die Werte in Photopeakposition an den gefundenen Indizes durch den Mittelwert
    Photopeakposition[indices] = average_value

    # Entferne die Werte in Feldstaerke und Photopeakposition an den gefundenen Indizes
    Feldstaerke = np.delete(Feldstaerke, indices)
    Photopeakposition = np.delete(Photopeakposition, indices)

    return Photopeakposition, Feldstaerke