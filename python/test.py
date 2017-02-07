import numpy as np
import matplotlib.pyplot as plt

import slr

sliceThickness = 5.
z = np.linspace(-3*sliceThickness,3*sliceThickness,400)

durationSLR = 3.1
rfSLR = slr.slr("se",400,3.55,durationSLR,flipAngle=160,filterType="ls")
rfSLR.GenerateRF()
rfUnscaledSLR = rfSLR.GetRF()
rfScaledSLR = rfSLR.GetRFScaled()
mxySLR = rfSLR.Simulate(sliceThickness,z)
timeSLR = np.linspace(0,durationSLR,400)


durationMsinc = 1.2
rfMsinc = slr.slr("smalltip",400,1.2,durationMsinc,flipAngle=160.)
rfMsinc.GenerateRF()
rfUnscaledMsinc = rfMsinc.GetRF()
rfScaledMsinc = rfMsinc.GetRFScaled()
mxyMsinc = rfMsinc.Simulate(sliceThickness,z,simulationType="se")
timeMsinc = np.linspace(0,durationMsinc,400)

plt.figure()
p1, = plt.plot(timeMsinc,rfScaledMsinc.real)
p2, = plt.plot(timeSLR,rfScaledSLR.real)
plt.legend([p1,p2],["m-sinc","SLR"])

plt.figure()
p1, = plt.plot(z,abs(mxyMsinc))
p2, = plt.plot(z,abs(mxySLR))
plt.legend([p1,p2],["m-sinc","SLR"])
plt.show()
