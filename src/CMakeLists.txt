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
alloc.hpp
cell.hpp
constOxyCell.hpp
constOxyTissue.hpp
coupler.hpp
createAndReadFiles.hpp
evalModel.hpp
model.hpp
morris.hpp
oxyCell.hpp
oxyTissue.hpp
rootSimulator.hpp
sensAn.hpp
simpcell.hpp
simulator.hpp
sobol.hpp
tcp.hpp
tissue.hpp
treatment.hpp
var.hpp
)

# RTsim basic source files
set(RTSIM_SOURCE_FILES
absOxyCell.cpp
alloc.cpp
cell.cpp
constOxyCell.cpp
constOxyTissue.cpp
coupler.cpp
createAndReadFiles.cpp
evalModel.cpp
model.cpp
morris.cpp
oxyCell.cpp
oxyTissue.cpp
rootSimulator.cpp
sensAnR.cpp
simulator.cpp
sobol.cpp
tcp.cpp
tissue.cpp
treatment.cpp
var.cpp
)

set(TOOLS
main.cpp
)

ADD_EXECUTABLE(RTsim ${RTSIM_SOURCE_FILES} ${TOOLS})



