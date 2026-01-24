cmake_minimum_required(VERSION 3.15)
project(velesquant LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare)
endif()

# Find Dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(QuantLib REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
# Handle pybind11 matching the Python version
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _pybind11_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmake_dir}")
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 QUIET)

# If Eigen3 not found via config, try standard include paths or environment variables
if(NOT EIGEN3_INCLUDE_DIR AND NOT TARGET Eigen3::Eigen)
    find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
        HINTS ${EIGEN3_ROOT} $ENV{EIGEN3_ROOT} /opt/homebrew/include/eigen3 /opt/homebrew/include
        PATH_SUFFIXES include/eigen3 include
    )
endif()

message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")

# OpenMP Check (set variables here, apply to targets later)
find_package(OpenMP)
set(USE_OPENMP_FALLBACK FALSE)
if(NOT OpenMP_FOUND AND APPLE)
    # Fallback for Homebrew OpenMP on macOS
    set(OMP_HOMEBREW_ROOT "")
    if(EXISTS "/opt/homebrew/opt/libomp")
        set(OMP_HOMEBREW_ROOT "/opt/homebrew/opt/libomp")
    elseif(EXISTS "/usr/local/opt/libomp")
        set(OMP_HOMEBREW_ROOT "/usr/local/opt/libomp")
    endif()

    if(OMP_HOMEBREW_ROOT)
        message(STATUS "Found Homebrew OpenMP at ${OMP_HOMEBREW_ROOT}")
        set(USE_OPENMP_FALLBACK TRUE)
        set(OMP_FALLBACK_INCLUDE "${OMP_HOMEBREW_ROOT}/include")
        set(OMP_FALLBACK_LIB "${OMP_HOMEBREW_ROOT}/lib/libomp.dylib")
    else()
         message(WARNING "OpenMP not found. Compilation may fail if code relies on <omp.h>")
    endif()
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

if(TARGET Eigen3::Eigen)
     target_link_libraries(velesquant_core PUBLIC Eigen3::Eigen)
endif()

# Apply OpenMP settings
if(OpenMP_FOUND)
    target_link_libraries(velesquant_core PRIVATE OpenMP::OpenMP_CXX)
elseif(USE_OPENMP_FALLBACK)
    target_compile_options(velesquant_core PUBLIC -Xpreprocessor -fopenmp)
    target_include_directories(velesquant_core PUBLIC "${OMP_FALLBACK_INCLUDE}")
    target_link_libraries(velesquant_core PUBLIC "${OMP_FALLBACK_LIB}")
endif()

# Create the Python Extension Module (modular bindings)
pybind11_add_module(_core
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

target_link_libraries(_core PRIVATE velesquant_core)
set_target_properties(_core PROPERTIES 
    INTERPROCEDURAL_OPTIMIZATION FALSE
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)
target_compile_options(_core PRIVATE -fno-lto -fvisibility=default)
target_link_options(_core PRIVATE -fno-lto)

# Install the module to the python package
install(TARGETS _core DESTINATION velesquant)

# Tests
add_subdirectory(tests_cpp)
