#
#@BEGIN LICENSE
#
# myscf by Psi4 Developer, a plugin to:
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

#
# Plugin Makefile generated by Psi4.
#
# You shouldn't need to modify anything in this file.
#

# The name of your plugin. Taken from the directory name.
NAME = $(shell basename `pwd`)

# C++ source files for your plugin. By default we grab all *.cc files.
CXXSRC = $(notdir $(wildcard *.cc))

# Flags that were used to compile Psi4.
CXX = /usr/bin/g++
CXXDEFS = -DFC_SYMBOL=2 -DHAVE_SYSTEM_NATIVE_LAPACK -DHAVE_SYSTEM_NATIVE_BLAS -DHAS_CXX11_VARIADIC_TEMPLATES -DHAS_CXX11_STATIC_ASSERT -DHAS_CXX11_SIZEOF_MEMBER -DHAS_CXX11_RVALUE_REFERENCES -DHAS_CXX11_LIB_REGEX -DHAS_CXX11_NULLPTR -DHAS_CXX11_LONG_LONG -DHAS_CXX11_LAMBDA -DHAS_CXX11_INITIALIZER_LIST -DHAS_CXX11_DECLTYPE -DHAS_CXX11_CSTDINT_H -DHAS_CXX11_CONSTEXPR -DHAS_CXX11_AUTO_RET_TYPE -DHAS_CXX11_AUTO -DHAS_CXX11_FUNC -DHAS_CXX11 -DVAR_MFDS -DSYS_DARWIN
CXXFLAGS = -DRESTRICT=__restrict__ -fPIC -std=c++11 -O3 -DNDEBUG -Wno-unused
INCLUDES = -I/Users/jeffschriber/Source/psi4/obj/src/lib -I/Users/jeffschriber/Source/psi4/src/lib -I/Users/jeffschriber/Source/psi4/include -I/Users/jeffschriber/Source/psi4/obj/include -I/Users/jeffschriber/Source/psi4/obj/boost/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/include
OBJDIR = /Users/jeffschriber/Source/psi4/obj

# Used to determine linking flags.
UNAME = $(shell uname)

# Need to link against Psi4 plugin library
PSIPLUGIN = -L$(OBJDIR)/lib -lplugin

DEPENDINCLUDE = $(notdir $(wildcard *.h*))

PSITARGET = $(NAME).so

# Start the compilation rules
default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -shared
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
    CXXFLAGS += -fno-common
endif

# The object files
BINOBJ = $(CXXSRC:%.cc=%.o)

%.o: %.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<

$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS) $(PSIPLUGIN)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d *.pyc *.test output.dat psi.timer.dat

