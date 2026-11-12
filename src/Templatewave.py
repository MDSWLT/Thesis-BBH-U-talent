import pycbc.waveform

class Templatewave:
    def __init__(self, m1, m2):
        self.mass1 = m1
        self.mass2 = m2
        self.hp, self.hc = pycbc.waveform.get_td_waveform(approximant='IMRPhenomD', mass1=m1, mass2=m2, distance=300, f_lower=20, delta_t=1 / 4096)
        self.hc = None