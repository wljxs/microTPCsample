# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(microTPCComsol)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
set(Eigen3_DIR "/home/liugx/Software/eigen-3.4.0/share/eigen3/cmake")
find_package(Eigen3 REQUIRED NO_MODULE)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR "Eigen3 not found!")
endif()

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(ElectronDict SHARED)

root_generate_dictionary(G__Electron
  Electron.h
  LINKDEF ElectronLinkDef.h
  MODULE ElectronDict)

target_link_libraries(ElectronDict PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree)
# ########################################################################################

# ---Build executable------------------------------------------------------------

add_executable(microTPC microTPC.C)
target_link_libraries(microTPC  Garfield::Garfield ElectronDict Eigen3::Eigen ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree)

add_executable(microTPCsimulation microTPCsimulation.C)
target_link_libraries(microTPCsimulation  Garfield::Garfield ElectronDict Eigen3::Eigen ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree)


# ---Copy all files locally to the build directory-------------------------------
# set(Detector microTPC)

# foreach(_file mplist.txt field_dr_300_d_3_dlc_450.txt xlabel.txt ylabel.txt mesh.mphtxt ar_95_iC4H10_5_1-200000.gas)
#   configure_file(${Detector}/${_file} ${_file} COPYONLY)
# endforeach()
