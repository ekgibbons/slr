"""Framework for generating, visualizing, and simulating SLR RF pulses

Instantiate appropriate class with filename.  Returned object acts like a
dictionary, with key-value pairs for each piece of metadata.
    import slr
    rf = slr.GenerateRF(numberPoints,
"""

__author__ = "Translated from C code from John Pauly by Eric Gibbons"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 01/30/2017 $"
__copyright__ = "(c) 2017 Eric Gibbons"
__license__ = "Python"

import matplotlib.pyplot as plt
import numpy as np

import slr_c

class slr:
    """ Generates the SLR pulse
    """
    
    def __init__(self, rf_type, numberPoints, tbw, duration=None,
                 inRipple=0.01, outRipple=0.01):
        """
        
        Arguments:
        - `rf_type`: 'ex', 'smalltip', 'se', 'sat','inv'
        - `numberPoints`: number of points in rf pulse
        - `tbw`: the time bandwidth
        - `duration`: the pulse duration [ms]
        - `inRipple`: the passband ripple (default is 0.001)
        - `outRipple`: the stopband ripple (default is 0.001)
        """

        availType = ["ex","smalltip","se","sat","inv"]

        try:
            pulse_index = availType.index(rf_type)
        except ValueError:
            print "ERROR:  please enter valid rf pulse type"
            print "\t'%s' is not valid" % rf_type
            print "\tOptions are: 'ex','smalltip','se','sat','inv'"

        self.rf_type = pulse_index
        self.numberPoints = numberPoints
        self.tbw = tbw
        self.duration = duration
        self.inRipple = inRipple
        self.outRipple = outRipple

        self.rf = slr_c.GenerateRF(self.numberPoints,self.tbw,self.rf_type,
                                   self.inRipple,self.outRipple)


    def RFScaleG(self,duration=None):
        """
        
        Arguments:
        - `time`: the pulse duration if not previously specified [ms]
        """

        if duration is not None:
            self.duration = duration


        self.rfScaled = self.rf*len(self.rf)/(self.duration*2*np.pi)
        self.rfScaled /= (4258./1.e3)

        return self.rfScaled
        
    def PlotRF(self,scaling=True):
        """
        
        Arguments:
        - `scaling`:
        """

        if scaling:
            self.RFScaleG()
            rfPlot = self.rfScaled
        else:
            rfPlot = self.rf

        plt.figure()
        plt.plot(rfPlot)
        plt.show()        
        
    def GetRF(self):
        """
        
        Arguments:

        """

        return self.rf
        
        
        
