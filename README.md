# slr - a Python implemenation of SLR pulse design
Python implementation of SLR RF design.  Much of this is based on the Pauly/LeRoux paper and the Pauly RF tools Matlab package. 

## Installation

This package does require some external packages, particularly the latest version of boost (v1.63).  At the moment this is installed by source.  Information can be found [here](http://www.boost.org/doc/libs/1_63_0/more/getting_started/unix-variants.html).

To use, you will first need to build the program.  Assuming you are on a Linux machine
```shell
cd slr
mkdir obj # build objects will be placed here
make
make install
```
The `install` step determines where you place the install path.  The install location is established by 
```shell
INSTALL_PATH = /location/of/python/path
```
## Usage
At the moment is is very rough.  An example script `test.py` is included to show a very simple pulse.  At the moment, the module has only one usable function to generate an RF pulse.  More to come later...


For a test case, see the `python/test.py` script for a quick example.  There is an SLR object slowly being built.  The documentation can be found by launching python in the command line and then
```python
>>>import slr
>>>help(slr)
```
to see that documentation.


The excitation seems to be working well, the other pulses...not as much.  Hopefully can get those working soon.
