/* SLR WRAPPER
 * 
 * This file is a wrapper for much of the C library.  This is written in 
 * C++ and is best compiled with the supplied Makefile.  You will also 
 * need the Boost Python (and Boost Python Numpy) libraries to have this 
 * work.  This generates the python module "rf_tools" that runs the same 
 * as any python module, though it is faster since the backend was written
 * and compiled in C.
 * 
 * Written by:  Eric Gibbons (Stanford University)
 * Written on:  2017.02.03
 * 
 */ 

#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "numpy_wrappers.hpp"

extern "C" {
#include "slr_tools.h"
}

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

#if _MSC_VER
using boost::uint8_t;
#endif

namespace py = boost::python;
namespace np = boost::python::numpy;

/* Wrapper for the Beta2Alpha function.  This function takes the beta polynomial
 * as a numpy array and calculates the alpha polynomial and returns it to python
 * as a complex numpy array.  Note:  the I/O is in double complex!
 */
np::ndarray SimulateRF(const np::ndarray &rfNumpy, const np::ndarray &gradientNumpy,
		       const np::ndarray &xNumpy, const std::string &rfType)
{
    int numberRF = (int)rfNumpy.shape(0);
    int numberGradient = (int)gradientNumpy.shape(0);
    int numberX = (int)xNumpy.shape(0);

    ASSERT_THROW((numberRF-numberGradient)==0,"RF and gradient lengths must match");

    double *rf = static_cast<double *>(malloc(2*numberRF*sizeof(double)));
    double *gradient = static_cast<double *>(malloc(numberGradient*sizeof(double)));
    double *x = static_cast<double *>(malloc(numberX*sizeof(double)));

    double *alpha = static_cast<double *>(calloc(2*(size_t)numberX,sizeof(double)));
    double *beta = static_cast<double *>(calloc(2*(size_t)numberX,sizeof(double)));
    double *profile = static_cast<double *>(malloc(2*numberX*sizeof(double)));

    std::string dtype;
    dtype = py::extract<std::string>(py::str(rfNumpy.get_dtype()));

    if (dtype.compare("complex128") == 0)
    {
    	CNumpy2CArray(rfNumpy, rf, numberRF);
    }
    else if (dtype.compare("float64") == 0)
    {
    	Numpy2CArray(rfNumpy, rf, numberRF);
    }

    Numpy2Array(gradientNumpy, gradient, numberGradient);
    Numpy2Array(xNumpy, x, numberX);

    abrx(rf,gradient,numberRF,x,numberX,alpha,beta);
    
    if (rfType.compare("ex") == 0)
    {
	ab2ex(profile,alpha,beta,numberX);
    }
    else if (rfType.compare("se") == 0)
    {
	ab2se(profile,beta,numberX);
    }
    else if (rfType.compare("inv") == 0)
    {
	ab2inv(profile,beta,numberX);
    }
    else if (rfType.compare("sat") == 0)
    {
	ab2inv(profile,beta,numberX);
    }
    else
    {
	throw std::runtime_error("Must choose proper simulation type!");
    }

    np::ndarray result = CArray2CNumpy(profile,numberX);

    free(rf);
    free(gradient);
    free(x);
    
    free(alpha);
    free(beta);
    free(profile);

    
    return result;
}


/* Wrapper for the Beta2Alpha function.  This function takes the beta polynomial
 * as a numpy array and calculates the alpha polynomial and returns it to python
 * as a complex numpy array.  Note:  the I/O is in double complex!
 */
np::ndarray Beta2AlphaWrap(const np::ndarray &betaNumpy)
{
    int numberPoints = (int)betaNumpy.shape(0);
    
    double *beta = static_cast<double *>(malloc(2*numberPoints*sizeof(double)));
    double *alpha = static_cast<double *>(malloc(2*numberPoints*sizeof(double)));

    std::string dtype = py::extract<std::string>(py::str(betaNumpy.get_dtype()));
    
    if (dtype.compare("complex128") == 0)
    {
	CNumpy2CArray(betaNumpy, beta, numberPoints);
    }
    else if (dtype.compare("float64") == 0)
    {
	Numpy2CArray(betaNumpy, beta, numberPoints);
    }

    Beta2Alpha(alpha, beta, numberPoints);
    
    np::ndarray result = CArray2CNumpy(alpha,numberPoints);

    free(beta);
    free(alpha);

    return result;
}


/* Wrapper for the Beta2RF function.  This function takes the beta polynomial as a 
 * numpy array and calculates the alpha polynomial.  It then takes beta and polynomial 
 * and applies the inverse SLR transform to find the RF, which is returned as to 
 * python as a numpy array.  Note:  the I/O is in double complex!
 */
np::ndarray Beta2RF(const np::ndarray &betaNumpy)
{
    int numberPoints = (int)betaNumpy.shape(0);
    
    double *beta = static_cast<double *>(malloc(2*numberPoints*sizeof(double)));
    double *alpha = static_cast<double *>(malloc(2*numberPoints*sizeof(double)));
    double *rf = static_cast<double *>(malloc(2*numberPoints*sizeof(double)));

    std::string dtype = py::extract<std::string>(py::str(betaNumpy.get_dtype()));

    if (dtype.compare("complex128") == 0)
    {
	CNumpy2CArray(betaNumpy, beta, numberPoints);
    }
    else if (dtype.compare("float64") == 0)
    {
	Numpy2CArray(betaNumpy, beta, numberPoints);
    }

    Beta2Alpha(alpha, beta, numberPoints);
    InverseSLR(rf, alpha, beta, numberPoints);
    
    np::ndarray result = CArray2CNumpy(rf,numberPoints);

    free(beta);
    free(alpha);
    free(rf);
    
    return result;
}

BOOST_PYTHON_MODULE(rf_tools)
{
    Py_Initialize();
    np::initialize();
    py::def("Beta2Alpha",Beta2AlphaWrap);
    py::def("Beta2RF",Beta2RF);
    py::def("SimulateRF",SimulateRF);
}
