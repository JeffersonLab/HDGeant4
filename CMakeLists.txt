cmake_minimum_required(VERSION 3.16)
project(HDGeant4 VERSION 1.0 LANGUAGES CXX)

# --- Standard C++ Detection (Matches your GCC logic) ---
set(CMAKE_CXX_STANDARD 11)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 11)
        set(CMAKE_CXX_STANDARD 20)
    elseif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
        set(CMAKE_CXX_STANDARD 17)
    elseif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 5)
        set(CMAKE_CXX_STANDARD 14)
    endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Find Dependencies ---
find_package(Geant4 REQUIRED ui_all vis_all gdml)
include(${Geant4_USE_FILE})

find_package(HDDM REQUIRED)

if(${Geant4_VERSION_MAJOR} LESS 10)
    set(G4_MAJOR_PADDED ${Geant4_VERSION_MAJOR})
    set(G4_MINOR_PADDED ${Geant4_VERSION_MINOR})
    set(G4_PATCH_PADDED ${Geant4_VERSION_PATCH})
else()
    set(G4_MAJOR_PADDED ${Geant4_VERSION_MAJOR})
    if(${Geant4_VERSION_MINOR} LESS 10)
        set(G4_MINOR_PADDED 0${Geant4_VERSION_MINOR})
    else()
        set(G4_MINOR_PADDED ${Geant4_VERSION_MINOR})
    endif()
    if(${Geant4_VERSION_PATCH} LESS 10)
        set(G4_PATCH_PADDED 0${Geant4_VERSION_PATCH})
    else()
        set(G4_PATCH_PADDED ${Geant4_VERSION_PATCH})
    endif()
endif()
set(G4FIXES_DIR_NAME "G4.${G4_MAJOR_PADDED}.${G4_MINOR_PADDED}.p${G4_PATCH_PADDED}fixes")
set(G4FIXES_DIR_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/${G4FIXES_DIR_NAME}")
if(IS_DIRECTORY ${G4FIXES_DIR_PATH})
    message(STATUS "Using G4fixes from ${G4FIXES_DIR_PATH}")
else()
    message(FATAL_ERROR "Cannot find G4fixes directory ${G4FIXES_DIR_PATH}, please create and try again.")
endif()

find_package(OpenGL REQUIRED)
find_path(MY_GLU_INCLUDE_DIR GL/glu.h 
    PATHS /usr/include /usr/local/include
    NO_DEFAULT_PATH
)
find_library(MY_GLU_LIBRARY NAMES GLU 
    PATHS /usr/lib64 /usr/lib
    NO_DEFAULT_PATH
)
message(STATUS "GLU Header: ${MY_GLU_INCLUDE_DIR}")
message(STATUS "GLU Library: ${MY_GLU_LIBRARY}")

set(OLD_LOG_LEVEL ${CMAKE_MESSAGE_LOG_LEVEL})
set(CMAKE_MESSAGE_LOG_LEVEL ERROR)
find_package(ROOT REQUIRED COMPONENTS Geom TMVA TreePlayer)
include(${ROOT_USE_FILE})
set(CMAKE_MESSAGE_LOG_LEVEL ${OLD_LOG_LEVEL})

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

find_package(Boost REQUIRED COMPONENTS python)

# For XERCES (required by HDDS and Geant4 GDML)
find_package(XercesC REQUIRED)

find_package(Diracxx REQUIRED)

find_package(JANA REQUIRED)

find_package(OpenMP REQUIRED)

# Handle Manual/Environment-based Dependencies (DANA/HALLD)
set(HALLD_RECON_HOME $ENV{HALLD_RECON_HOME})
set(HDDS_HOME $ENV{HDDS_HOME})
set(JANA_HOME $ENV{JANA_HOME})

# Options to be selected at build time
option(USE_DIRACXX "Use Dirac++ library" ON)
if(USE_DIRACXX)
    add_definitions(-DUSING_DIRACXX)
endif()

option(BYPASS_DRAWING_CLIPPED_VOLUMES "Bypass code that clips hidden volumes in visualization" ON)
if(BYPASS_DRAWING_CLIPPED_VOLUMES)
    add_definitions(-DLAYERED_GEOMETRY_PICKING_EXTENSIONS)
endif()

set(REDUCE_OPTIMIZATION_OF_CDC "1" CACHE STRING "CDC optimization reduction level")
add_definitions(-DREDUCE_OPTIMIZATION_OF_CDC=${REDUCE_OPTIMIZATION_OF_CDC})

option(G4VIS_BUILD_OPENGL_DRIVER "Build the OpenGL driver" ON)
if(G4VIS_BUILD_OPENGL_DRIVER)
    add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)
endif()

option(G4VIS_BUILD_OPENGLX_DRIVER "Build the OpenGLX driver" ON)
if(G4VIS_BUILD_OPENGLX_DRIVER)
    add_definitions(-DG4VIS_BUILD_OPENGLX_DRIVER)
endif()

option(G4MULTITHREADED "Build a multithreaded executable" ON)
if(G4MULTITHREADED)
    add_definitions(-DG4MULTITHREADED)
endif()

option(USE_G4VIS "Enable Geant4 visualization" ON)
if(USE_G4VIS)
    add_definitions(-DG4VIS_USE)
endif()

option(USE_G4UI "Enable the Geant4 interactive user interface" ON)
if(USE_G4UI)
    add_definitions(-DG4UI_USE)
endif()

# Conditional Macro for DIRC
if(EXISTS "${HALLD_RECON_HOME}/src/libraries/DIRC/DDIRCPmtHit.cc")
    add_definitions(-DDIRCTRUTHEXTRA)
endif()

# --- Include Directories ---
include_directories(
    src
    src/${G4FIXES_DIR_NAME}
    src/G4debug
    ${HDDS_HOME}
    ${HALLD_RECON_HOME}/$ENV{BMS_OSNAME}/include
    ${JANA_HOME}/include
    ${XercesC_INCLUDE_DIRS}
)

# --- Define Libraries ---

# 1. HDDS Library
add_library(hdds SHARED 
    ${HDDS_HOME}/XString.cpp 
    ${HDDS_HOME}/XParsers.cpp 
    ${HDDS_HOME}/hddsCommon.cpp
)

# 2. G4fixes Library
file(GLOB G4FIXES_SRC "src/${G4FIXES_DIR_NAME}/*.cc" "src/G4debug/*.cc")
add_library(G4fixes SHARED src/G4fixes.cc ${G4FIXES_SRC})
target_compile_definitions(G4fixes PRIVATE ${Geant4_DEFINITIONS})
target_link_libraries(G4fixes ${Geant4_LIBRARIES} Boost::python Python3::Python ${MY_GLU_LIBRARY})

# 3. Cobrems Library
file(GLOB COBREMS_SRC "src/Cobrems*.cc")
add_library(cobrems SHARED ${COBREMS_SRC})
target_link_libraries(cobrems ${Geant4_LIBRARIES} Boost::python Python3::Python)

# 4. Main hdgeant4 Library
option(WITH_G4VERBOSE "Enable Geant4 verbose output" OFF)
option(WITH_GEOM_DEBUG "Enable heavy geometry debugging" OFF)
option(WITH_TRACKING_DEBUG "Enable heavy tracking debugging" OFF)
if(WITH_G4VERBOSE)
    target_compile_definitions(hdgeant4 PUBLIC
       G4VERBOSE
    )
endif()
if(WITH_GEOM_DEBUG)
    target_compile_definitions(hdgeant4 PRIVATE
       G4GEOMETRY_DEBUG
       G4CSGDEBUG
       G4DIVDEBUG
    )
endif()
if(WITH_TRACKING_DEBUG)
    target_compile_definitions(hdgeant4 PRIVATE
       G4_TRACKING_DEBUG
       G4DEBUG_FIELD
       G4DEBUG_NAVIGATION
       G4DEBUG_PARALLEL_WORLD_PROCESS
       G4_USE_G4BESTUNIT_FOR_VERBOSE
       G4DEBUG_TRANSPORT
       G4DEBUG_VIS_OGL
    )
endif()

file(GLOB HDGEANT4_SRC "src/*.cc")
file(GLOB JANA2_LIBS "${JANA_HOME}/lib/lib*.so")
list(REMOVE_ITEM HDGEANT4_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/CobremsGeneration.cc" "${CMAKE_CURRENT_SOURCE_DIR}/src/G4fixes.cc")
add_library(hdgeant4 SHARED ${HDGEANT4_SRC})
target_compile_definitions(hdgeant4 PUBLIC ${Geant4_DEFINITIONS})
target_compile_definitions(hdgeant4 PUBLIC
    G4UI_USE
    G4VIS_USE
    G4UI_USE_EXECUTIVE
    G4VIS_USE_EXECUTIVE
    G4VIS_BUILD_OPENGL_DRIVER
    G4VIS_BUILD_OPENGLX_DRIVER
)
target_link_directories(hdgeant4 PUBLIC 
    "${HALLD_RECON_HOME}/$ENV{BMS_OSNAME}/lib"
    "$ENV{EVIOROOT}/lib" 
    "$ENV{SQLITECPP_HOME}/lib"
    "$ENV{SQLITE_HOME}/lib"
    "$ENV{CCDB_HOME}/lib"
)
target_link_libraries(hdgeant4 PUBLIC
    hdds G4fixes cobrems
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES}
    ${XercesC_LIBRARIES}
    Diracxx::Dirac
    JANA::jana2_shared_lib
    BCAL CDC DANA DIRC ECAL FCAL FDC FMWPC CCAL ECAL
    RF TAGGER PAIR_SPECTROMETER START_COUNTER TOF TRD
    ANALYSIS HDGEOMETRY PID TAC TPOL TRIGGER CERE
    TRACKING KINFITTER TTAB DAQ EVENTSTORE
    HDDM
    HDDM::xstream
    evioxx evio
    ccdb mysqlclient
    SQLiteCpp sqlite3
    GLU
)

# --- Executable ---
add_executable(hdgeant4_bin src/hdgeant4.cc) # Assuming main is here
set_target_properties(hdgeant4_bin PROPERTIES OUTPUT_NAME hdgeant4)
target_link_libraries(hdgeant4_bin hdgeant4)

# --- Utils ---
add_subdirectory(src/utils)

install(TARGETS hdgeant4_bin
        RUNTIME DESTINATION bin)

install(TARGETS hdgeant4 cobrems G4fixes
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(TARGETS hdgeant4 DESTINATION g4py/HDGeant4)
install(TARGETS cobrems  DESTINATION g4py/Cobrems)
install(TARGETS G4fixes  DESTINATION g4py/G4fixes)

install(FILES 
        g4py/cobrems.py 
        g4py/hdgeant4.py 
        DESTINATION g4py)
install(FILES g4py/Cobrems/__init__.py DESTINATION g4py/Cobrems)
install(FILES g4py/G4fixes/__init__.py DESTINATION g4py/G4fixes)
install(FILES g4py/HDGeant4/__init__.py DESTINATION g4py/HDGeant4)
