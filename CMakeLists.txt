cmake_minimum_required(VERSION 3.15)
project(velesquant LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Dependencies
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/quantlib")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/boost")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(QuantLib REQUIRED)
find_package(Boost REQUIRED)
find_package(pybind11 REQUIRED)
find_package(OpenMP)

if(OpenMP_FOUND)
    message(STATUS "OpenMP found")
else()
    message(STATUS "OpenMP not found only with find_package, trying homebrew path")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor" "-fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    include_directories("/opt/homebrew/opt/libomp/include")
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Include Directories
include_directories(
    include
    ${QuantLib_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Collect Source Files
file(GLOB_RECURSE SOURCES
    "src/local_vol/*.cpp"
    "src/pde_solvers/*.cpp"
    "src/models/*.cpp"
)
list(FILTER SOURCES EXCLUDE REGEX "test_simplex.cpp")

# Create the Core C++ Library
# We exclude the original main/interface files from the glob if they existed (checked earlier, seemingly not relevant for lib)
add_library(velesquant_core STATIC ${SOURCES})

target_link_libraries(velesquant_core
    PRIVATE
    ${QuantLib_LIBRARIES}
    ${Boost_LIBRARIES}
)
if(OpenMP_FOUND)
    target_link_libraries(velesquant_core PRIVATE OpenMP::OpenMP_CXX)
else()
    # Fallback for manual OpenMP setup if find_package failed
    if(OpenMP_omp_LIBRARY)
        target_link_libraries(velesquant_core PRIVATE "${OpenMP_omp_LIBRARY}")
    endif()
endif()
if(OpenMP_CXX_FLAGS)
    target_compile_options(velesquant_core PRIVATE ${OpenMP_CXX_FLAGS})
endif()

# Create the Python Extension Module
pybind11_add_module(_core
    bindings/bindings.cpp
)

target_link_libraries(_core PRIVATE velesquant_core)

# Install the module to the python package
install(TARGETS _core DESTINATION velesquant)

# Tests
add_subdirectory(tests_cpp)
