cmake_minimum_required(VERSION 3.15)
project(velesquant LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare)
endif()
# Find Dependencies
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/quantlib")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/boost")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/eigen")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(QuantLib REQUIRED)
find_package(Boost REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 QUIET)
# If find_package didn't set EIGEN3_INCLUDE_DIR, use fallback
if(NOT EIGEN3_INCLUDE_DIR OR EIGEN3_INCLUDE_DIR STREQUAL "")
    # Fallback for Homebrew Eigen
    if(APPLE)
        set(EIGEN3_INCLUDE_DIR "/opt/homebrew/opt/eigen/include/eigen3")
    endif()
endif()
message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")
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
# Include Directories moved to target setup
# include_directories() removed

# Collect Source Files
file(GLOB_RECURSE SOURCES
    "src/volatility/*.cpp"
    "src/pde_solvers/*.cpp"
    "src/models/*.cpp"
    "src/engines/*.cpp"
    "src/instruments/*.cpp"
    "src/numerics/*.cpp"
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
target_include_directories(velesquant_core PUBLIC include)
target_include_directories(velesquant_core SYSTEM PUBLIC
    ${QuantLib_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
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

# Create the Python Extension Module (modular bindings)
pybind11_add_module(native
    bindings/module.cpp
    bindings/bind_enums.cpp
    bindings/bind_exceptions.cpp
    bindings/bind_utility.cpp
    bindings/models/bind_heston.cpp
    bindings/models/bind_hull_white.cpp
    bindings/models/bind_sabr.cpp
    bindings/models/bind_local_vol.cpp
    bindings/models/bind_swaption.cpp
    bindings/models/bind_cms.cpp
    bindings/models/bind_basket.cpp
    bindings/models/bind_hybrid.cpp
    bindings/models/bind_trees.cpp
    bindings/pde_solvers/bind_pde_solvers.cpp
)

target_link_libraries(native PRIVATE velesquant_core)

# Install the module to the python package
install(TARGETS native DESTINATION velesquant)

# Tests
add_subdirectory(tests_cpp)
