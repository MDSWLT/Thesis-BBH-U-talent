import numpy as np
import matplotlib.pyplot as plt
import pycbc.waveform
from matplotlib import pyplot as plt
import math


class Templatewave:
    def __init__(self, m1, m2, d):
        self.mass1 = m1
        self.mass2 = m2
        self.distance = d
        self.hp , self.hc = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4P', mass1=m1, mass2=m2, disance=d, f_lower=20, delta_t=1/4096)

#Generen van een Templatebank
templatebank = []

for m1 in range(20):
    for m2 in range(20):
        templatewave = Templatewave(11+(m1/10), 11+(m2/10), 300)
        templatebank.append(templatewave)
        print("wave gegenereerd")

hp1, hc1 = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4P', mass1=12, mass2=12, disance=300, f_lower=20, delta_t=1/4096)
hp2, hc2 = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4P', mass1=10.5, mass2=13.5, disance=300, f_lower=20, delta_t=1/4096)
hp3, hc3 = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4P', mass1=14, mass2=10, disance=300, f_lower=20, delta_t=1/4096)

array = np.zeros((40,40))

#Wavematching met een voorlopig fautieve formule
for t in templatebank:
    if np.shape(hp1.numpy()) == np.shape(t.hp.numpy()):
        overlap = (np.sum(hp1.numpy()*t.hp.numpy())/math.sqrt(np.sum(hp1.numpy()**2)*np.sum(t.hp.numpy()**2)))
        print(f"Een massa van {t.mass1} en een massa van {t.mass2} geven een overlap van {overlap}")
        array[round((t.mass1-11)*10), round((t.mass2-10)*10)] = overlap
heatmap = plt.imshow(array)

fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)
ax[0].plot(np.arange(len(hp1))/4096 - 27, hp1)
ax[1].plot(np.arange(len(hp2))/4096 - 27, hp2)
ax[2].plot(np.arange(len(hp3))/4096 - 27, hp3)
plt.show()
