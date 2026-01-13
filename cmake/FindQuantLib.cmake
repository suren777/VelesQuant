# FindQuantLib.cmake

find_path(QuantLib_INCLUDE_DIR ql/quantlib.hpp
    HINTS
    /opt/homebrew/include
    /opt/homebrew/opt/quantlib/include
    /usr/local/include
    ${QuantLib_DIR}/include
)

find_library(QuantLib_LIBRARY NAMES QuantLib
    HINTS
    /opt/homebrew/lib
    /opt/homebrew/opt/quantlib/lib
    /usr/local/lib
    ${QuantLib_DIR}/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QuantLib
    REQUIRED_VARS QuantLib_LIBRARY QuantLib_INCLUDE_DIR
)

if(QuantLib_FOUND)
    set(QuantLib_LIBRARIES ${QuantLib_LIBRARY})
    set(QuantLib_INCLUDE_DIRS ${QuantLib_INCLUDE_DIR})
    
    if(NOT TARGET QuantLib::QuantLib)
        add_library(QuantLib::QuantLib UNKNOWN IMPORTED)
        set_target_properties(QuantLib::QuantLib PROPERTIES
            IMPORTED_LOCATION "${QuantLib_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${QuantLib_INCLUDE_DIR}"
        )
    endif()
endif()
