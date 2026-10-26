import numpy as np
import pycbc.waveform
from matplotlib import pyplot as plt
import pycbc.noise
import pycbc.psd.analytical
import math

# aan te roepen functie waarin de onbekende massa bepaald kan worden, aan de hand van de timeseries van de hp
def wavematch(hp, iterations):
    class Templatewave:
        def __init__(self, m1, m2, d):
            self.mass1 = m1
            self.mass2 = m2
            self.distance = d
            self.hp , self.hc = pycbc.waveform.get_td_waveform(approximant='IMRPhenomD', mass1=m1, mass2=m2, distance=d, f_lower=20, delta_t=1/4096)

    mstart = 4
    overlaps = []
    masses = []
    # ntemplates moet voor de eerste iteratie 400 zijn om veel massa's te checken, maar dit is bij de volgende iteraties niet nodig dus wordt deze kleiner.
    ntemplates = 400
    for i in range(iterations):
        # Het genereren van de templatebank
        templatebank = []
        for m in range(ntemplates):
            templatewave = Templatewave(mstart + (m / 10**(i+1)), 12, 300)
            if np.shape(hp1.numpy()) == np.shape(templatewave.hp.numpy()):
                templatebank.append(templatewave)
        # Het matchen zelf
        for t in templatebank:
            overlap = (np.sum(hp.numpy() * t.hp.numpy()) / math.sqrt(np.sum(hp.numpy() ** 2) * np.sum(t.hp.numpy() ** 2)))
            overlaps.append(overlap)
            masses.append(t.mass1)
        mclosest = masses[overlaps.index(max(overlaps))]
        mstart = mclosest - 1/10**(i+1)
        ntemplates = 20
    return masses, overlaps

hp1, hc1 = pycbc.waveform.get_td_waveform(approximant='IMRPhenomD', mass1=26.690, mass2=12, distance=440, f_lower=20, delta_t=1/4096)

deltaf = 1/(len(hp1)/4096)
length = int(4096/2/deltaf +1)
psd = pycbc.psd.analytical.AdvVirgo(length, deltaf, 20)
noise = pycbc.noise.noise_from_psd(length=len(hp1), delta_t=1/4096, psd=psd, seed=100)

signaal = np.add(hp1, noise)

masses, overlaps = wavematch(signaal, 3)
plt.plot(masses, overlaps, 'o')
plt.show()

# SNR berekenen
def bereken_SNR(no_noise, noise):
    SNR_squared = np.sum(no_noise.numpy()**2) / np.sum(noise.numpy()**2)
    return np.sqrt(SNR_squared)

# aanroepen SNR
SNR = bereken_SNR(hp1, noise)
print(SNR)
