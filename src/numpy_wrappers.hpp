/* numpy_wrappers.hpp
 * 
 * Just some basic wrapping functions to move numpy
 * formatted array to standard C/C++ arrays
 *
 */


#ifndef __NUMPY_WRAPPERS_HPP__
#define  __NUMPY_WRAPPERS_HPP__

namespace py = boost::python;
namespace np = boost::python::numpy;

np::ndarray Array2Numpy(double *array, int length);
np::ndarray CArray2CNumpy(double *array, int length);
void Numpy2Array(const np::ndarray &arrayNumpy,
		 double *arrayOut,
		 int length);
void Numpy2CArray(const np::ndarray &arrayNumpy,
		  double *arrayOut,
		      int length);
void CNumpy2CArray(const np::ndarray &arrayNumpy,
		  double *arrayOut,
		  int length);

#endif
