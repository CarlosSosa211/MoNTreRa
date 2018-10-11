# General CMake configuration
################################################################################

cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE OFF)
add_compile_options(-std=c++11)
################################################################################
# General Project configuration
################################################################################

# project name
project(RTsim)

################################################################################
# Source code declarations
################################################################################

# RTsimulator basic header files

set(RTSIM_HEADER_FILES
absOxyCell.hpp
cell.hpp
constOxyCell.hpp
constOxyTissue.hpp
coupler.hpp
model.hpp
oxyCell.hpp
oxyTissue.hpp
rootSimulator.hpp
simpcell.hpp
simulator.hpp
tissue.hpp
treatment.hpp
)

# RTsim basic source files
set(RTSIM_SOURCE_FILES
absOxyCell.cpp
cell.cpp
constOxyCell.cpp
constOxyTissue.cpp
coupler.cpp
model.cpp
oxyCell.cpp
oxyTissue.cpp
rootSimulator.cpp
simulator.cpp
tissue.cpp
treatment.cpp
)

set(TOOLS
main.cpp
)

ADD_EXECUTABLE(RTsim ${RTSIM_SOURCE_FILES} ${TOOLS})



