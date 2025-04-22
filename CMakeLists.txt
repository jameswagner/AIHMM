cmake_minimum_required(VERSION 3.10)
project(AIHMM VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add options for build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

# Find required packages - specifically for Ubuntu
find_package(Boost REQUIRED COMPONENTS program_options)
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost not found. Please install Boost:\n   sudo apt-get install libboost-program-options-dev")
endif()

# Source files in new structure
set(SOURCES
    # Core files
    src/core/Chromosome.cpp
    src/core/ExpressionInfo.cpp
    src/core/Individual.cpp
    src/core/SNP.cpp
    
    # HMM files
    src/hmm/HMMSNP.cpp
    src/hmm/constants.c
    
    # I/O files
    src/io/FileReader.cpp
    
    # Main file
    src/main.cpp
)

# Include directories
include_directories(
    ${Boost_INCLUDE_DIRS}
    src/core
    src/hmm
    src/io
)

# Create executable
add_executable(AIHMM ${SOURCES})

# Link libraries
target_link_libraries(AIHMM PRIVATE
    ${Boost_LIBRARIES}
)

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Installation rules
install(TARGETS AIHMM DESTINATION bin)
install(FILES 
        models/ceu8StateHMMParams.txt 
        DESTINATION share/AIHMM)

# Copy sample parameter file to build directory
file(COPY models/ceu8StateHMMParams.txt DESTINATION ${CMAKE_BINARY_DIR})

# Print configuration summary
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Boost Include Dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "  Boost Libraries: ${Boost_LIBRARIES}")
message(STATUS "  Output Directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}") 