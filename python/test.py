import numpy as np
import matplotlib.pyplot as plt

import slr
import slr_c

rf = slr.slr('ex',200,4.,duration=2.)
rf_pulse = rf.GetRF() # return the actual pulse values
rf.PlotRF()

