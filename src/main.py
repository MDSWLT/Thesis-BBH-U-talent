import numpy as np
import pycbc.waveform
from matplotlib import pyplot as plt
import pycbc.noise
import pycbc.psd.analytical
import math
import statistics

import pickle
def load(path):
    with open(path, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                return

import pathlib
datapad = f"{pathlib.Path(__file__).parent.parent.absolute()}/basistemplatebank"

from Templatewave import Templatewave

# De eerste 400 templatewaves worden meteen aan het begin gegenereerd, zodat dit rekenwerk maar een keer hoeft te gebeuren
basistemplatebank = []
for m in range(700):
    basistemplatebank.append(Templatewave(4 + m/10, 22))

# Aan te roepen functie voor het genereren van 2n+1 templatewaves met massa's in de spreiding mass±n/10^decimal (decimal heeft
# de waarde van de positie van het te veranderen cijfer achter de komma. Bijvoorbeeld: 0 -> 1.0, 2.0, 3.0 etc. en 1 -> 0.1, 0.2, 0.3 etc.)
def maketemplates(n, mass, decimal):
    templatebank = []
    mstart = mass - (n/(10**decimal))
    for i in range(2*n + 1):
        templatebank.append(Templatewave(mstart + (i/(10**decimal)), 22))
    return templatebank

# Aan te roepen functie voor het matchen van een golfsignaal in de vorm van een timeseries met templatewaves (in de vorm
# van lijst met class instances van Templatewave). Er wordt gekeken naar de hp. De output is een lijst met alle testmassa's
# en een lijst met alle overlapwaardes
def wavematch(hp, templatebank):
    masses = []
    overlaps = []
    for t in templatebank:
        if np.shape(hp.numpy()) == np.shape(t.hp.numpy()):
            overlap = abs(np.sum(hp.numpy() * t.hp.numpy()) / math.sqrt(np.sum(hp.numpy() ** 2) * np.sum(t.hp.numpy() ** 2)))
            masses.append(t.mass1)
            overlaps.append(overlap)
    return masses, overlaps

# Aan te roepen functie die van een hp timeseries en een waarde voor de hoeveelheid getallen achter de komma
# rond de beste match een lijst met testmassa's en een lijst met bijbehorende overlaps maakt.
def findmass(hp, decimal):
    masses, overlaps = wavematch(hp, basistemplatebank)
    bestmass = masses[overlaps.index(max(overlaps))]
    for i in range(decimal):
        templatebank = maketemplates(9, bestmass, i+1)
        extramasses, extraoverlaps = wavematch(hp, templatebank)
        bestmass = extramasses[extraoverlaps.index(max(extraoverlaps))]
        for m in extramasses:
            masses.append(m)
        for o in extraoverlaps:
            overlaps.append(o)
    return masses, overlaps

for d in range(10):
    hp1, hc1 = pycbc.waveform.get_td_waveform(approximant='IMRPhenomD', mass1=58, mass2=22, distance=100*(d+1), f_lower=20, delta_t=1/4096)
    deltaf = 1/(len(hp1)/4096)
    length = int(4096/2/deltaf +1)
    psd = pycbc.psd.analytical.AdvVirgo(length, deltaf, 20)

    gevonden_massas = []
    for i in range(20):
        noise = pycbc.noise.noise_from_psd(length=len(hp1), delta_t=1/4096, psd=psd, seed=i)
        signaal = np.add(hp1, noise)
        masses, overlaps = findmass(signaal, 3)
        bestmass = masses[overlaps.index(max(overlaps))]
        gevonden_massas.append(bestmass)
    print(f"Afstand: {100*(d+1)} Mpc")
    print(f"De gemiddelde gevonden massa is: {statistics.mean(gevonden_massas)} met een standaarddeviatie: {statistics.stdev(gevonden_massas)}")

# SNR berekenen
def bereken_SNR(no_noise, noise):
    SNR_squared = np.sum(no_noise.numpy()**2) / np.sum(noise.numpy()**2)
    return np.sqrt(SNR_squared)

# aanroepen SNR
#SNR = bereken_SNR(hp1, noise)
#print(SNR)
