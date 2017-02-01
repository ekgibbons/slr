#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>


extern "C" {
#include "slr_design.h"
}


#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

#if _MSC_VER
using boost::uint8_t;
#endif

namespace py = boost::python;
namespace np = boost::python::numpy;

np::ndarray Array2Numpy(float *array, int length)
{
    py::tuple shape = py::make_tuple(length);
    np::dtype dt1 = np::dtype::get_builtin<std::complex<double> >(); 
    np::ndarray result = np::zeros(shape, dt1);

    for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
    {
	result[ii] = (double)array[ii];
    }
    
    return result;
}

np::ndarray GenerateRF(const int &nSamples,
		       const float &sliceThick,
		       const float &duration,
		       const float &amp,
		       const float &inRipple,
		       const float &outRipple,
		       const float &type)
{
    float *rf = static_cast<float *>(malloc(nSamples*sizeof(float)));

    gen_slr_rf(rf,nSamples,sliceThick,duration,amp,inRipple,
	       outRipple,type);
    
    np::ndarray result = Array2Numpy(rf,nSamples);
    
    return result;
}



BOOST_PYTHON_MODULE(slr_c)
{
    Py_Initialize();
    np::initialize();
    py::def("GenerateRF",GenerateRF);

}
