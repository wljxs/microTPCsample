cmake_minimum_required(VERSION 3.14)
project(BeamAnalysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------
# 查找 ROOT
# ------------------------------
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
    message(FATAL_ERROR "ROOT not found!")
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR "Eigen3 not found!")
endif()

set(SOURCES
    src/main.cpp
    src/Event/DetectorFrame.cpp
    src/Event/EventDisplay.cpp
    src/Detector/Detector.cpp
    src/Detector/Planar.cpp
    src/Detector/Cylinder.cpp
    src/Detector/DetectorFactory.cpp
    src/Algorithm/AlgorithmFactory.cpp
    src/Algorithm/algorithms/WaveformProcessor.cpp
    src/Algorithm/algorithms/ClusterBuilder.cpp
    src/Algorithm/algorithms/ClusterReconstructor.cpp
    src/Algorithm/AnalysisUtils.cpp
    src/Algorithm/WaveformAnalysis.cpp
    src/Script/Base/RawDataParser.cpp
    src/Script/Base/ScriptFactory.cpp
    src/Script/Base/ScriptManager.cpp
    src/Script/ExampleScript.cpp
    src/Script/TrackAnalysisScript.cpp
    src/Script/DeconvAnalysisScript.cpp
    src/Script/DUTAnalysisScript.cpp
    src/Script/BadStripScript.cpp
    src/Script/ChannelNoiseScript.cpp
    src/Script/RawWaveformDumpScript.cpp
)

add_executable(BeamAnalysis ${SOURCES})

target_include_directories(BeamAnalysis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} 
)

target_link_libraries(BeamAnalysis
    ${ROOT_LIBRARIES}
    stdc++fs
)

# ------------------------------
# 共享库（供 ROOT 宏 .x 加载使用，含反卷积、波形分析等功能）
# ------------------------------
# 从 SOURCES 中移除 main.cpp，保留其余所有编译单元
set(LIB_SOURCES ${SOURCES})
list(REMOVE_ITEM LIB_SOURCES src/main.cpp)

add_library(BeamAnalysisLib SHARED ${LIB_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/BeamAnalysis_dict.cxx)

target_include_directories(BeamAnalysisLib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(BeamAnalysisLib
    ${ROOT_LIBRARIES}
    stdc++fs
)

# 头文件和 LinkDef
set(DI_HDR ${CMAKE_CURRENT_SOURCE_DIR}/include/Event/DataModel.h)
set(LINKDEF  ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)

ROOT_GENERATE_DICTIONARY(BeamAnalysis_dict ${CMAKE_CURRENT_SOURCE_DIR}/include/Event/DataModel.h 
                        MODULE BeamAnalysis 
                        LINKDEF ${LINKDEF})

# ------------------------------
# 仿真分析可执行程序
# ------------------------------
add_executable(SimulationAnalysis src/simulationanalysis.C)

target_include_directories(SimulationAnalysis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(SimulationAnalysis
    BeamAnalysisLib
    ${ROOT_LIBRARIES}
)

# ------------------------------
# Tikhonov 正则化参数扫描
# ------------------------------
add_executable(LambdaScan src/lambda_scan.C)

target_include_directories(LambdaScan PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(LambdaScan
    BeamAnalysisLib
    ${ROOT_LIBRARIES}
)

add_executable(test src/test.C)

target_include_directories(test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Detector
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Event
    ${CMAKE_CURRENT_SOURCE_DIR}/include/EventDisplay
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Script
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(test
    BeamAnalysisLib
    ${ROOT_LIBRARIES}
)