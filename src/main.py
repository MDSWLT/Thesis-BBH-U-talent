eerste_massa = 26
tweede_massa = 26


import numpy as np
import pycbc.waveform
import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import pycbc.noise
import pycbc.psd.analytical
import math
import statistics
import pandas as pd

import pickle
def load(path):
    with open(path, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                return

import pathlib
datapad = pathlib.Path(__file__).parent.parent.absolute()

from Templatewave import Templatewave

totaltemplatebank = []
usedmasses = []

# Aan te roepen functie voor het genereren van templatewaves met massa's in de spreiding mass(1 en 2)±n/10^decimal (decimal heeft
# de waarde van de positie van het te veranderen cijfer achter de komma. Bijvoorbeeld: 0 -> 1.0, 2.0, 3.0 etc. en 1 -> 0.1, 0.2, 0.3 etc.)
def maketemplates(n, mass1, mass2, decimal):
    mstart1 = mass1 - (n/(10**decimal))
    mstart2 = mass2 - (n/(10**decimal))
    templatebank = []
    for i in range(2*n + 1):
        for j in range(2*n + 1):
            m1 = mstart1 + (i/(10**decimal))
            m2 = mstart2 + (j/(10**decimal))
            if f'{m1}{m2}' in usedmasses:
                templatebank.append(totaltemplatebank[usedmasses.index(f'{m1}{m2}')])
            else:
                wave = Templatewave(m1, m2)
                templatebank.append(wave)
                totaltemplatebank.append(wave)
                usedmasses.append(f'{m1}{m2}')
                print(f'Nieuwe golf gegenereerd ({m1}, {m2})')
    return templatebank

# Aan te roepen functie voor het matchen van een golfsignaal in de vorm van een timeseries met templatewaves (in de vorm
# van lijst met class instances van Templatewave). Er wordt gekeken naar de hp. De output is twee lijsten met alle testmassa's
# en een lijst met alle overlapwaardes
def wavematch(hp, matchtemplatebank):
    masses1 = []
    masses2 = []
    overlaps = []
    for t in matchtemplatebank:
        if np.shape(hp.numpy()) == np.shape(t.hp.numpy()):
            overlap = abs(np.sum(hp.numpy() * t.hp.numpy()) / math.sqrt(np.sum(hp.numpy() ** 2) * np.sum(t.hp.numpy() ** 2)))
            masses1.append(t.mass1)
            masses2.append(t.mass2)
            overlaps.append(overlap)
    return masses1, masses2, overlaps

# Aan te roepen functie die van een hp timeseries en een waarde voor de hoeveelheid getallen achter de komma
# (minimaal 1) rond de beste match een lijst met testmassa's en een lijst met bijbehorende overlaps maakt.
def findmass(hp, decimal):
    decimal -= 1
    masses1, masses2, overlaps = wavematch(hp, list(load(f'{datapad}/templatebank.pickle')))
    bestmass1 = masses1[overlaps.index(max(overlaps))]
    bestmass2 = masses2[overlaps.index(max(overlaps))]
    for i in range(decimal):
        templatebank = maketemplates(6, bestmass1, bestmass2, i+2)
        extramasses1, extramasses2, extraoverlaps = wavematch(hp, templatebank)
        bestmass1 = extramasses1[extraoverlaps.index(max(extraoverlaps))]
        bestmass2 = extramasses2[extraoverlaps.index(max(extraoverlaps))]
        for m in extramasses1:
            masses1.append(m)
        for m in extramasses2:
            masses2.append(m)
        for o in extraoverlaps:
            overlaps.append(o)
    return masses1, masses2, overlaps

with pd.ExcelWriter(f"{datapad}/resultaten({eerste_massa} en {tweede_massa}).xlsx", engine="openpyxl") as writer:
    gemiddelde_massas1 = []
    gemiddelde_massas2 = []
    snrs = []
    standaarddeviaties1 = []
    standaarddeviaties2 = []

    for d in range(40):
        distance = 8000/(d+1)
        hp, hc = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4P', mass1=eerste_massa, mass2=tweede_massa, distance=distance, f_lower=20, delta_t=1/4096)
        deltaf = 1/(len(hp)/4096)
        length = int(4096/2/deltaf +1)
        psd = pycbc.psd.analytical.AdvVirgo(length, deltaf, 20)

        gevonden_overlaps = []
        gevonden_massas1 = []
        gevonden_massas2 = []
        for i in range(75):
            noise = pycbc.noise.noise_from_psd(length=len(hp), delta_t=1/4096, psd=psd)
            signaal = np.add(hp, noise)
            masses1, masses2, overlaps = findmass(signaal, 3)
            bestoverlap = max(overlaps)
            bestmass1 = masses1[overlaps.index(bestoverlap)]
            bestmass2 = masses2[overlaps.index(bestoverlap)]
            gevonden_overlaps.append(bestoverlap)
            gevonden_massas1.append(bestmass1)
            gevonden_massas2.append(bestmass2)
            print(bestmass1, bestmass2)

        tussenresultaten = pd.DataFrame({
            "Gevonden massa 1": gevonden_massas1,
            "Gevonden massa 2": gevonden_massas2,
            "Overlap": gevonden_overlaps,
        })

        tussenresultaten.to_excel(writer, sheet_name=f"d={distance}", index=False)

        snr = pycbc.filter.sigma(hp, psd, low_frequency_cutoff=20, high_frequency_cutoff=2048)
        snrs.append(snr)

        gemiddelde_massa1 = statistics.mean(gevonden_massas1)
        gemiddelde_massas1.append(gemiddelde_massa1)

        standaarddeviatie_massa1 = statistics.stdev(gevonden_massas1)
        standaarddeviaties1.append(standaarddeviatie_massa1)

        gemiddelde_massa2 = statistics.mean(gevonden_massas2)
        gemiddelde_massas2.append(gemiddelde_massa2)

        standaarddeviatie_massa2 = statistics.stdev(gevonden_massas2)
        standaarddeviaties2.append(standaarddeviatie_massa2)

        print(f"SNR: {snr}\n"
            f"M1: {gemiddelde_massa1} met een standaarddeviatie: {standaarddeviatie_massa1}\n"
            f"M2: {gemiddelde_massa2} met een standaarddeviatie: {standaarddeviatie_massa2}")

    eindresultaten = pd.DataFrame({
        "SNR": snrs,
        "Massa 1": gemiddelde_massas1,
        "Standaarddeviatie 1": standaarddeviaties1,
        "Massa 2": gemiddelde_massas2,
        "Standaarddeviatie 2": standaarddeviaties2,
    })

    eindresultaten.to_excel(writer, sheet_name=f"Eindresultaten", index=False)
