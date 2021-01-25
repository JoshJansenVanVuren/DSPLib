cmake_minimum_required(VERSION 3.0.0)
project(DSPLib VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_library(FFT lib/fft.cpp lib/fft.h)
add_library(STFT lib/stft.cpp lib/stft.h)
add_library(Doppler lib/doppler.cpp lib/doppler.h)
add_library(Filter lib/filter.cpp lib/filter.h)
add_executable(RunRadar src/main.cpp)

INCLUDE_DIRECTORIES(lib/ )
target_link_libraries(RunRadar PRIVATE FFT STFT Filter Doppler)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include(CPack)
