/* numpy_wrappers.cpp
 * 
 * Just some basic wrapping functions to move numpy
 * formatted array to standard C/C++ arrays
 *
 */


#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <math.h>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "numpy_wrappers.hpp"

#if _MSC_VER
using boost::uint8_t;
#endif

namespace py = boost::python;
namespace np = boost::python::numpy;

np::ndarray Array2Numpy(double *array, int length)
{
    py::tuple shape = py::make_tuple(length);
    np::dtype dt1 = np::dtype::get_builtin<double>(); 
    np::ndarray result = np::zeros(shape, dt1);

    for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
    {
	result[ii] = array[ii];
    }
    
    return result;
}


np::ndarray CArray2CNumpy(double *array, int length)
{
    std::complex<double> j = std::complex<double>(0,1.);
    
    py::tuple shape = py::make_tuple(length);
    np::dtype dt1 = np::dtype::get_builtin<std::complex<double> >(); 
    np::ndarray result = np::zeros(shape, dt1);

    for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
    {
	result[ii] = array[2*ii] + j*array[2*ii + 1];
    }
    
    return result;
}

void Numpy2Array(const np::ndarray &arrayNumpy,
		 double *arrayOut,
		 int length)
{

    std::vector<double> tempVector(length);
    
    std::copy(reinterpret_cast<double*>(arrayNumpy.get_data()), 
	      reinterpret_cast<double*>(arrayNumpy.get_data())+length,
	      tempVector.begin());

    for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
    {
	arrayOut[ii] = (double)tempVector[ii];
    }
}

void Numpy2CArray(const np::ndarray &arrayNumpy,
		  double *arrayOut,
		  int length)
{

    std::vector<double> tempVector(length);
    
    std::copy(reinterpret_cast<double*>(arrayNumpy.get_data()), 
	      reinterpret_cast<double*>(arrayNumpy.get_data())+length,
	      tempVector.begin());

    for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
    {
	arrayOut[2*ii] = (double)tempVector[ii];
    }
}


void CNumpy2CArray(const np::ndarray &arrayNumpy,
		  double *arrayOut,
		  int length)
{

    std::vector<std::complex<double> > tempVector(length);
    
    std::copy(reinterpret_cast<std::complex<double>* >(arrayNumpy.get_data()), 
	      reinterpret_cast<std::complex<double>* >(arrayNumpy.get_data())+length,
	      tempVector.begin());

    for (unsigned int ii = 0; ii < (unsigned int) length; ii++)
    {
	printf("%f + j%f\n",tempVector[ii].real(),tempVector[ii].imag());
	arrayOut[2*ii] = (double)tempVector[ii].real();
	arrayOut[2*ii+1] = (double)tempVector[ii].imag();
    }
}
