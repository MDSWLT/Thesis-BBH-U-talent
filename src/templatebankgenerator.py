# LET OP!! Het runnen van dit programma duurt lang! (~10 uur op een relatief sterke computer)

import pickle
import pathlib

from Templatewave import Templatewave

path = pathlib.Path(__file__).parent.parent.absolute()

with open(f"{path}/templatebank.pickle", "wb") as f:
    for m in range(301):
        print(round(6 + m / 10, 1))
        for n in range(301-m):
            pickle.dump(Templatewave(round(6 + m/10, 1), round(36 - (n/10), 1)), f)
