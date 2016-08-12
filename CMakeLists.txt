# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10)
project (decimatedDistributions)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

# If ROOTSYS is defined in the usual way then this should be sufficient to find the ROOT CMake things
if(DEFINED $ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# If not then we can access the variables the ROOT CMake stuff directly, fingers crossed
else()
  execute_process(COMMAND root-config --etcdir OUTPUT_VARIABLE ROOT_ETC_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ROOT_ETC_DIR}/cmake)
endif()

if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
elseif(DEFINED ENV{ARA_UTIL_INSTALL_DIR})
    message("ARA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ARA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ARA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ARA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()

find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit Minuit2)
message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
add_definitions(${ROOT_CXX_FLAGS})
link_directories(${UTIL_LIB_DIR})

set(ANITA_LIBS AnitaEvent AnitaAnalysis AnitaAnalysisTools AnitaCorrelator RootFftwWrapper)

#add_library(CalibRuns SHARED CalibRuns.cxx)
#target_link_libraries(CalibRuns ${ROOT_LIBRARIES} ${ANITA_LIBS})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wshadow -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Werror")

set(BINARIES reconstruct)

FOREACH(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable(${binary} ${binary}.cxx)
  target_link_libraries(${binary} ${ZLIB_LIBRARIES} ${ANITA_LIBS} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})  
ENDFOREACH(binary)



