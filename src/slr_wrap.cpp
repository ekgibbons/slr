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
    np::dtype dt1 = np::dtype::get_builtin<double>(); 
    np::ndarray result = np::zeros(shape, dt1);

    for (unsigned int ii = 0; ii < (unsigned int)length; ii++)
    {
	result[ii] = (double)array[ii];
    }
    
    return result;
}

np::ndarray GenerateRF(const int &nSamples,
		       const float &tbw,
		       const int &ptype,
		       const float &inRipple,
		       const float &outRipple)
{
    float *rf = static_cast<float *>(malloc(nSamples*sizeof(float)));

    gen_slr_rf(rf,nSamples,tbw,ptype,inRipple,outRipple);
    
    np::ndarray result = Array2Numpy(rf,nSamples);
    
    return result;
}



BOOST_PYTHON_MODULE(slr_c)
{
    Py_Initialize();
    np::initialize();
    py::def("GenerateRF",GenerateRF);

}
