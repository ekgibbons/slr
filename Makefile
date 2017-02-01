TARGET = slr_c
# OBJECTS = slr_wrap.o slr_design.o remez.o
# #  = slr_design.h remez.h

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES   = $(wildcard $(SRCDIR)/*.c) #$(wildcard $(SRCDIR)/*.cpp)
SOURCES_CPP = $(wildcard $(SRCDIR)/*.cpp)
INCUDES   = $(wildcard $(SRCDIR)/*.h) #$(wildcard $(SRCDIR)/*.hpp)
INCUDES_HPP = $(wildcard $(SRCDIR)/*.hpp)
OBJECTS_C    = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
# OBJECTS_C    += $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
OBJECTS_CPP  = $(SOURCES_CPP:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# location of the Python header files
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/local/lib

# boost libraries
LIBRARY_PATH = -L$(BOOST_LIB)
PYTHON_LIBS = -lpython$(PYTHON_VERSION) -lboost_python -lboost_numpy 

LD_LIBS = $(PYTHON_LIBS)


CPPFLAGS = -O3 -Wall -Wconversion -fPIC 
CXXFLAGS = -std=c++11 -I$(PYTHON_INCLUDE) -I$(BOOST_INC)
# CFLAGS = 
LDFLAGS = -shared -Wl,--export-dynamic 
LDFLAGS += -Wl,--unresolved-symbols=report-all


$(TARGET): $(OBJECTS_C) $(OBJECTS_CPP)
	$(CXX) $(LDFLAGS) $(OBJECTS_C) $(OBJECTS_CPP) $(LIBRARY_PATH) $(LD_LIBS) -o $(TARGET).so

$(OBJECTS_C): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CPPFLAGS) -c -o $@ $< 

$(OBJECTS_CPP): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

install:
	cp $(TARGET).so ~/python/modules/mri

clean:
	rm -rf *~ *.o *.so
	rm -rf $(OBJECTS_C)
	rm -rf $(OBJECTS_CPP)
	rm -rf $(TARGET).so

