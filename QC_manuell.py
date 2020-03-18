# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:57:43 2020

@author: bm2384
"""
# Import Section
import time
import numpy as np
from datetime import datetime
from datetime import timedelta
import glob
import os
import soundfile as sf
import matplotlib.pyplot as plt
# Set all parameters for import
t = time.time()
path = 'X:/Meereskunde/Unterwasserschall/BSH-Projekt-UBA-PIMO/Messstelle_FeBe/FeBe_WAV_201901-201906/FeBe_WAV_Slot_A/Data/'
dateDeployment = datetime.strptime('2019-01-29-14-13-00', '%Y-%m-%d-%H-%M-%S')
dateRecovery = datetime.strptime('2019-06-11-13-13-00', '%Y-%m-%d-%H-%M-%S')
DutyCycleDauer = timedelta(minutes=15)
DutyCycleAbstand = timedelta(minutes=60)
SollSampleRate = 32000

# Scan durch den angegebenen Ordner
os.chdir(path)
NoRealData = path + '/NoRealData/'
Samplerates = []
FileLength = []
ProzClippedSamples = []
Mittelwerte = []
Q_1 = []
Q_2 = []
Q_3 = []
Q_4 = []
Q_5 = []
Q_6 = []
Q_7 = []
Q_8 = []

try:
    os.mkdir(NoRealData)
except:
    print('Ordner NoRealData existiert schon')
f = open('QC_log.txt', 'w+')
f.write('Protokoll des Qualitätschecks \n \n'
        'Qualitätslegende:\n'
        '1: Gut \n'
        '2: Aussortiert - vor Deployment \n'
        '3: Aussortiert - nach Recovery \n'
        '4: Aussortiert - samplerate zu niedrig \n'
        '5: Aussortiert - Anteil geclippter Daten zu hoch \n'
        '6: Aussortiert - File ist 1% kürzer als Solllänge \n'
        '7: Aussortiert - File ist 1% länger als Solllänge \n'
        '8: Aussortiert - DC - Offset \n \n \n')


f.write('         filename                  | aussortiert? | Qualität \n'
        '--------------------------------------------------------------- \n')
FileListeAll = glob.glob("*.wav")
i = 0
for file in glob.glob("*.wav"):
    i = i + 1
    print('Bin bei file ' + str(i) + ' von ' + str(len(FileListeAll)))
    # Verschieben aller Daten, die davor osder danach aufgenommen wurden
    zeit = file[16:-4]
    zeit = datetime.strptime(zeit, '%Y%m%d_%H%M%S')
    signal, fs = sf.read(file)
    Samplerates.append(fs)
    FileLength.append(len(signal) / fs)  # Länge des files in Sekunden
    Mittelwerte.append(np.mean(signal))
    # Check for clipping
    ClippedSamples = np.where(abs(signal) == 1)
    AnteilClippedSamples = len(ClippedSamples) / len(signal) * 100
    ProzClippedSamples.append(AnteilClippedSamples)
    if zeit < dateDeployment.replace(second=0, minute=0) + DutyCycleAbstand:
        print('Noch vor Deployment')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    2     \n')
        Q_2.append(1)
    elif zeit > dateRecovery.replace(second=0, minute=0) - DutyCycleAbstand:
        print('Schon nach Recovery')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    3     \n')
        Q_3.append(1)
    elif fs < SollSampleRate:
        print('Samplerate hat sich geändert auf: ' + str(fs))
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    4     \n')
        Q_4.append(1)
    elif AnteilClippedSamples > 0.1:
        print(str(AnteilClippedSamples) + ' % clipped')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    5     \n')
        Q_5.append(1)
    elif len(signal)/fs < 60 * DutyCycleDauer.seconds/60 * 0.99:    # 1% zu kurz
        print('File ist zu kurz')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    6     \n')
        Q_6.append(1)
    elif len(signal)/fs > 60 * DutyCycleDauer.seconds/60 * 1.01:    # 1% zu lang
        print('File ist zu lang')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    7     \n')
        Q_7.append(1)
    elif abs(np.mean(signal)) > 0.01:
        print('DC- Offset ist größer als 1%')
        os.rename(path + file, NoRealData + file)
        f.write(file + '|      Yes     |    8     \n')
        Q_8.append(1)
    else:
        print('Alles gut')
        f.write(file + '|      No      |    1     \n')
        Q_1.append(1)
f.write('--------------------------------------------------------------- \n \n \n \n')       

# Herausfinden der tatsächlichen vorhandenen files
FileListe = glob.glob("*.wav")
AnzahlFilesIsso = len(FileListe)+1
print('Es sind insgesamt: ' + str(AnzahlFilesIsso) + ' files vorhanden')
f.write('Es sind insgesamt: ' + str(AnzahlFilesIsso) + ' files vorhanden \n')
# Zeiten = []
# for i in enumerate(FileListe):
#     Zeiten.append(datetime.strptime(FileListe[int(i[0])][16:-4], '%Y%m%d_%H%M%S'))
# Herausfinden der theoretisch nötigen files
Start = dateDeployment.replace(second=0, minute=0) + DutyCycleAbstand
Ende = dateRecovery.replace(second=0, minute=0) - DutyCycleAbstand
duration = Ende - Start
AnzahlFilesTheo = (duration.total_seconds()/60)/(DutyCycleAbstand.seconds/60)
print('Es sollten: ' + str(AnzahlFilesTheo) + ' vorhanden sein')
f.write('Es sollten: ' + str(AnzahlFilesTheo) + ' vorhanden sein \n')
Differenz_files = AnzahlFilesTheo - AnzahlFilesIsso
print('Die Differenz beträgt: ' + str(Differenz_files) + ' files')
f.write('Die Differenz beträgt: ' + str(Differenz_files) + ' files \n \n \n')

# Visualisierung der Samplerates - gerader Strich idealerweise
# plt.plot(Samplerates)
# plt.show()
# plt.title('Samplerates in allen .wav files')
# plt.xlabel('File Index')
# plt.ylabel('Samplerate')
# plt.grid(True)

# Zählen warum was aussortiert wurde
f.write(' Qualität | Anzahl Files \n')
f.write('------------------------ \n')
f.write('   1      |     ' + str(len(Q_1)) + '     \n')
f.write('   2      |     ' + str(len(Q_2)) + '     \n')
f.write('   3      |     ' + str(len(Q_3)) + '     \n')
f.write('   4      |     ' + str(len(Q_4)) + '     \n')
f.write('   5      |     ' + str(len(Q_5)) + '     \n')
f.write('   6      |     ' + str(len(Q_6)) + '     \n')
f.write('   7      |     ' + str(len(Q_7)) + '     \n')
f.write('   8      |     ' + str(len(Q_8)) + '     \n')

f.write('Insgesamt wurden ' + str(len(FileListeAll)) + ' gescannt \n')
elapsed = time.time() - t
f.write('Das alles hat ' + str(elapsed) + 'gedauert')
f.close()

# Duty Cycle ermitteln
# Check wann erstes und wann letztes file weggeschrieben wurden
# audioinfo und check länge check samplerate check clipping check sound floor 
# angabe wieviel %geclippt sind - rausgeschrieben muss <0.1%