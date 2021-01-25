# DSPLib
A relatively simple and raw C++ implementation of a few DSP routines. Completed as a self exercise. The library can filter and compute both STFT and FFTs (only in powers of 2), and do some basic bandpass filtering. The routines can be accessed through the classes, as documented in docs/...

## To compile
```
cd build
make
```

## Project Structure
```
├── bin                     # Compiled binary for example program
├── build                   # CMake and make files
├── docs                    # Doxygen documentation files
├── lib                     # Libraries for different DSP Routines
│   ├── doppler.cpp
│   ├── doppler.h
│   ├── fft.cpp
│   ├── fft.h
│   ├── filter.cpp
│   ├── filter.h
│   ├── stft.cpp
│   ├── stft.h
├── src                     # Contains an example run through of the library
├── test                    # Unit testing      
├── CMakeLists.txt          # CMake file for make file creation
└── README.md         
```

## Project Status
There are a few additions I would like to add to the project:
 * Unit-tests
 * Exception handling
 * Expand on handled types
 * Reading from a WAV file
 * etcetera

## Project Environment
 * c++ 98 or later
 * CMake 3.0.0
 * GNU Make 4.1