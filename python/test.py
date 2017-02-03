import numpy as np
import matplotlib.pyplot as plt

import slr

rf = slr.slr('ex',50,4.,1.2)
rf.GenerateRF()
rfScaled = rf.GetRFScaled()
rfUnscaled = rf.GetRF()

plt.figure()
plt.plot(rfScaled.real)
plt.plot(rfUnscaled.real)
plt.show()
