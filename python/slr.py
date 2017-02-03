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
from scipy import fftpack

import rf_tools
import fir


class slr:
    """ Generates the SLR pulse
    """
    
    def __init__(self, rfType, numberPoints, tbw, duration,
                 flipAngle=None,inRipple=0.01, outRipple=0.01):
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
            pulse_index = availType.index(rfType)
        except ValueError:
            print "ERROR:  please enter valid rf pulse type"
            print "\t'%s' is not valid" % rf_type
            print "\tOptions are: 'ex','smalltip','se','sat','inv'"

        self.rfType = rfType # pulse_index
        self.numberPoints = numberPoints
        self.tbw = tbw
        self.duration = duration
        self.flipAngle = flipAngle
        self.d1e = inRipple
        self.d2e = outRipple
        
        # self.rf = slr_c.GenerateRF(self.numberPoints,self.tbw,self.rf_type,
        #                            self.inRipple,self.outRipple)

    def __DInf(self,d1,d2):
        """Calculates dinf
        
        Arguments:
        - `d1`: passband ripple
        - `d2`: stopband ripple
        """

        a1 = 5.309e-3
        a2 = 7.114e-2
        a3 = -4.761e-1
        a4 = -2.66e-3
        a5 = -5.941e-1
        a6 = -4.278e-1
    
        l10d1 = np.log10(d1)
        l10d2 = np.log10(d2)
        
        di = ((a1*l10d1*l10d1+a2*l10d1+a3)*l10d2
                   + (a4*l10d1*l10d1+a5*l10d1+a6))

        return di

    def __RFScaleG(self):
        """
        
        Arguments:
        - `time`: the pulse duration if not previously specified [ms]
        """

        self.rfScaled = self.rf*len(self.rf)/(self.duration*2*np.pi)
        self.rfScaled /= (4258./1.e3)

    def GenerateRF(self):
        """Takes the input parameters from when the class was
           instantiated and generates the RF pulse through 
           SLR design
        
        Arguments:
        - `self`:
        """

        # determine pulse parameters based on what type of pulse it is
        if self.rfType is "ex":
            bsf = np.sqrt(0.5)
            d1 = np.sqrt(self.d1e/2)
            d2 = self.d2e/np.sqrt(2)
        elif self.rfType is "smalltip":
            bsf = 1
            d1 = self.d1e
            d2 = self.d2e
        elif self.rfType is "se":
            bsf = 1
            d1 = self.d1e/4
            d2 = np.sqrt(self.d2e)
        elif self.rfType is "sat":
            bsf = np.sqrt(0.5)
            d1 = self.d1e/2
            d2 = np.sqrt(self.d2e)
        elif self.rfType is "inv":
            bsf = 1
            d1 = self.d1e/8
            d2 = np.sqrt(0.5*self.d2e)
        else:
            print "ERROR:  please choose an appropriate pulse type"
            pass

        di = self.__DInf(d1,d2)

        # transition band width
        w = di/self.tbw

        # frequency bands
        f = np.array([0, (1-w)*(self.tbw/2), (1+w)*(self.tbw/2),
                      (float(self.numberPoints)/2)])/(float(self.numberPoints)/2)

        m = np.array([1,1,0,0])
        w = np.array([1,d1/d2])

        # build the beta polynomial
        betaPolynomial = fir.ls(self.numberPoints,f,m,w)

        # determine scaling factor
        if self.flipAngle is not None:
            # find betaProfile, to better approximate the ripples, zero-pad
            # the fft.  
            betaProfile = fftpack.fft(betaPolynomial,n=4*len(betaPolynomial))
            bsf = np.sin(flipAngle/2)/np.amax(abs(betaProfile))

        betaPolynomial *= bsf

        # insure that the betaPolynomial is complex...probably won't be...but
        # this is okay.
        betaPolynomial = betaPolynomial.astype(complex)

        self.rf = rf_tools.Beta2RF(betaPolynomial)
        self.__RFScaleG()
        
    def GetRF(self):
        """Returns the unscaled RF waveform

        Arguments:
        - `self`:
        """

        return self.rf        
        
    def GetRFScaled(self):
        """Returns the scaled RF waveform
        
        Arguments:
        - `self`:
        """

        return self.rfScaled
        
        # plt.figure()
        # plt.plot(betaPolynomial.real)
        # plt.plot(self.rf.real)
        # plt.show()





        
    # def PlotRF(self,scaling=True):
    #     """
        
    #     Arguments:
    #     - `scaling`:
    #     """

    #     if scaling:
    #         self.RFScaleG()
    #         rfPlot = self.rfScaled
    #     else:
    #         rfPlot = self.rf

    #     plt.figure()
    #     plt.plot(rfPlot)
    #     plt.show()        
        
    # def GetRF(self):
    #     """
        
    #     Arguments:

    #     """

    #     return self.rf
        
    # def __Beta2Alpha(self,beta):
    #     """Takes the beta polynomial and converts it to the alpha
    #        polynomial.  Usually done in C routine, but that isn't
    #        the algorithm isn't that complicated...
        
    #     Arguments:
    #     - `beta`:
    #     """

    #     betaSize = len(beta)
    #     padSize = 2**(int(np.log2(betaSize) - np.log2(betaSize)%1))*2*16
    #     betaFT = fftpack.fft(beta,n=padSize)

    #     betaFT /= np.sqrt(np.amax(abs(betaFT)))

    #     plt.figure()
    #     plt.plot(abs(betaFT))
    #     plt.show()
        
    #     alphaMag = np.sqrt(1 - abs(betaFT)**2)
    #     alphaFT = np.sqrt(np.log(alphaMag))

    #     alphaFT = fftpack.ifft(alphaFT,n=padSize)
    #     alphaFT[1:padSize/2-1] *= 2
    #     alphaFT[padSize/2+1:padSize] = 0
        
    #     alphaFT = fftpack.ifft(alphaFT,n=padSize)
    #     alphaFT /= padSize

    #     phase = -alphaFT.imag
    #     alphaFT = alphaMag*np.exp(1j*phase)

    #     alpha = fftpack.ifft(alphaFT,n=padSize)
    #     alpha = alpha[:betaSize]/padSize
    #     alpha = alpha[::-1]

    #     print alpha
        
        
