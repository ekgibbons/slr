import matplotlib.pyplot as plt

import slr_c


rf = slr_c.GenerateRF(100,0.5,4e-3,1,0.01,0.01,1)

print rf

plt.figure()
plt.plot(rf)
plt.show()
