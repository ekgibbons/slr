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

import rf_tools
import fir

class slr:
    """ Generates the SLR pulse
    """
    
    def __init__(self, rfType, numberPoints, tbw, duration,
                 flipAngle=None,filterType="pm",inRipple=0.01, outRipple=0.01):
        """
        
        Arguments:
        - `rf_type`: 'ex', 'smalltip', 'se', 'sat','inv'
        - `numberPoints`: number of points in rf pulse
        - `tbw`: the time bandwidth
        - `duration`: the pulse duration [ms]
        - `flipAngle`: the flip angle in degrees (default is chosen by 
                       pulse type if not specified here
        - `filterType`: the filter for the beta (default is 'pm')
                        -> 'pm': parks mcclellan
                        -> 'ls': least squares
        - `inRipple`: the passband ripple (default is 0.001)
        - `outRipple`: the stopband ripple (default is 0.001)
        """

        availType = ["ex","smalltip","se","sat","inv"]

        try:
            pulse_index = availType.index(rfType)
        except RuntimeError:
            print "ERROR:  please enter valid rf pulse type"
            print "\t'%s' is not valid" % rfType
            print "\tOptions are: 'ex','smalltip','se','sat','inv'"

        if rfType is "smalltip" and filterType is not "smalltip":
            filterType = "smalltip"
        
        
        if rfType is not "smalltip" and filterType is "smalltip":
            raise RuntimeError("smalltip filters not available for this pulse")

        
        self.rfType = rfType # pulse_index
        self.numberPoints = numberPoints
        self.tbw = tbw
        self.duration = duration
        self.flipAngle = flipAngle
        self.filterType = filterType
        self.d1e = inRipple
        self.d2e = outRipple

    def _DInf(self,d1,d2):
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

    def _RFScaleG(self):
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
            raise RuntimeError("please choose an appropriate pulse type")

        if self.flipAngle is not None:
            bsf = np.sin(self.flipAngle/2.*np.pi/180.)
            
        di = self._DInf(d1,d2)

        # transition band width
        w = di/self.tbw

        # frequency bands
        f = np.array([0, (1-w)*(self.tbw/2), (1+w)*(self.tbw/2),
                      (float(self.numberPoints)/2)])/(float(self.numberPoints)/2)/2
        w = np.array([1,d1/d2])

        # build the beta polynomial
        if self.filterType is "ls":
            m = np.array([1,1,0,0])
            betaPolynomial = fir.ls(self.numberPoints,f,m,w)
        elif self.filterType is "smalltip":
            if self.flipAngle is None:
                print "please specify flip angle for small tip"
                pass
            bsf = self.flipAngle/180.*np.pi
            betaPolynomial = fir.msinc(self.numberPoints,self.tbw)
        elif self.filterType is "pm":
            m = np.array([1,0])
            betaPolynomial = fir.remez(self.numberPoints,f,m,w)
        elif self.filterType is "min" or self.filterType is "max":
            n2 =2*self.numberPoints - 1
            di = 0.5*self._DInf(2*d1,0.5*d2*d2)
            w = di/self.tbw
            f = np.array([0, (1-w)*(self.tbw/2), (1+w)*(self.tbw/2),
                    (float(self.numberPoints)/2)])/(float(self.numberPoints)/2)/2
            w = np.array([1,2*d1/(0.5*d2*d2)])
            m = np.array([1.,0.])
            betaPolynomial = fir.remez(n2,f,m,w)
            betaPolynomial = fir.MinPhase(betaPolynomial)
            if self.filterType is "min":
                betaPolynomial = betaPolynomial[::-1]
        else:
            raise RuntimeError("please choose the appropriate filter type")
            
        betaPolynomial *= bsf

        # insure that the betaPolynomial is complex...probably won't be...but
        # this is okay.
        betaPolynomial = betaPolynomial.astype(complex)

        if self.filterType is "smalltip":
            self.rf = betaPolynomial
        else:
            self.rf = rf_tools.Beta2RF(betaPolynomial)
            
        self._RFScaleG()
        
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

    def Simulate(self,sliceThickness,z,simulationType=None):
        """
        
        Arguments:
        - `self`:
        - `pulseWidth`:
        - `x`:
        """

        if simulationType is None:
            if self.rfType is "smalltip":
                raise RuntimeError("simulation type needed")
            else:
                simulationType = self.rfType
            
        sliceThicknessCM = sliceThickness/10
        zCM = z/10
        dt = self.duration/self.numberPoints
        
        bandwidth = self.tbw/self.duration
        gradientAmplitude = bandwidth/sliceThicknessCM;
        gradient = 2*np.pi*gradientAmplitude*dt*np.ones((self.numberPoints))

        profile = rf_tools.SimulateRF(self.rf.astype(complex),
                                      gradient,zCM,simulationType)
        
        return profile
        
        
        
