import numpy as np
import matplotlib.pyplot as plt

import slr

rf = slr.slr('se',400,4.,4.,filterType='min')
rf.GenerateRF()
rfScaled = rf.GetRFScaled()
rfUnscaled = rf.GetRF()

plt.figure()
# plt.plot(rfScaled.real)
plt.plot(rfUnscaled.real)
plt.show()
