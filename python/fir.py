import numpy as np

def ls(numtaps, bands, desired, weight):
    """
    FIR filter design using least-squares error minimization.

    This calculates the coefficients for a Type 1 or II filter. The
    filter type is automatically determined based on whether or not 
    numtaps is even or odd.

    Calculate the filter coefficients for the linear-phase finite
    impulse response (FIR) filter which has the best approximation
    to the desired frequency response described by `bands` and
    `desired` in the least squares sense (i.e., the integral of the
    weighted mean-squared error within the specified bands is
    minimized).

    Design borrows heavily from the Matlab implementation (firls.m).

    Parameters
    ----------
    numtaps : int
        The number of taps in the FIR filter.  `numtaps` must be odd.
    bands : array_like
        A monotonic nondecreasing sequence containing the band edges in
        Hz. All elements must be non-negative and less than or equal to
        the Nyquist frequency given by `nyq`.
    desired : array_like
        A sequence the same size as `bands` containing the desired gain
        at the start and end point of each band.
    weight : array_like, optional
        A relative weighting to give to each band region when solving
        the least squares problem. `weight` has to be half the size of
        `bands`.

    Returns
    -------
    coeffs : ndarray
        Coefficients of the optimal (in a least squares sense) FIR filter.


    """

    N = np.copy(numtaps)
    F = np.copy(bands)
    M = np.copy(desired)
    W = np.copy(weight)

    Nodd = N % 2

    F /= 2
    W = np.sqrt(W)
    dF = np.diff(F)
    
    L = (N - 1)/2

    m = np.arange(0,L+1) + 0.5*(1 - Nodd)
    k = m
    
    I11, I22 = np.meshgrid(k,m,)

    I1 = I11 + I22
    I2 = - I11 + I22

    if Nodd:
        k = k[1:]
        b0 = 0

    b = np.zeros_like(k)
    G = np.zeros_like(I1)

    for s in np.arange(0,len(F),2):
        m = (M[s+1] - M[s])/(F[s+1] - F[s]) # slope
        b1 = M[s] - m*F[s]
        
        if Nodd:
            b0 += (b1*(F[s+1]-F[s]) + m/2*(F[s+1]**2 - F[s]**2)*abs(W[s/2])**2)
        
        b += ((m/(4*np.pi**2)*(np.cos(2*np.pi*k*F[s+1]) - np.cos(2*np.pi*k*F[s]))/(k**2))
              * abs(W[s/2])**2)
        
        b += ((F[s+1]*(m*F[s+1]+b1)*np.sinc(2*k*F[s+1])
               - F[s]*(m*F[s]+b1)*np.sinc(2*k*F[s]))
               * abs(W[s/2])**2)

        
        G += ((0.5*F[s+1]*(np.sinc(2*I1*F[s+1]) + np.sinc(2*I2*F[s+1]))
                 - 0.5*F[s]*(np.sinc(2*I1*F[s]) + np.sinc(2*I2*F[s])))
                * abs(W[s/2])**2)
        
    if Nodd:
        b = np.hstack((b0, b))

    
        
    a = np.linalg.pinv(G).dot(b)

    if Nodd:
        h = np.hstack((a[L:0:-1]/2, a[0], a[1:L+1]/2))
    else:
        h = 0.5*np.hstack((a[::-1],a))

    coeffs = h
        
    return coeffs

    
