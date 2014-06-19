#---------------------------------------------------
# OpenMMMBPol REFERENCE Platform
#
# Creates OpenMM library, base name=OpenMMMBPolReference.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMMBPolReference[_d].dll
#   OpenMMMBPolReference[_d].lib
#   OpenMMMBPolReference_static[_d].lib
# Unix:
#   libOpenMMMBPolReference[_d].so
#   libOpenMMMBPolReference_static[_d].a
#----------------------------------------------------

SET(OPENMM_REFERENCE_LIBRARY_NAME OpenMMMBPolReference)


# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_MBPOL_SOURCE_SUBDIRS .)

SET(SHARED_TARGET ${OPENMM_REFERENCE_LIBRARY_NAME})

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (MSVC)
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (MSVC)

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)


# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/  include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp                        ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)
SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/src)

# Create the library

INCLUDE_DIRECTORIES(${REFERENCE_INCLUDE_DIR})

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}               ${API_INCLUDE_FILES})

SET(OPENMM_LIBRARY_NAME OpenMM)

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(MAIN_OPENMM_LIB ${OPENMM_LIBRARY_NAME}_d)
ELSE (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(MAIN_OPENMM_LIB ${OPENMM_LIBRARY_NAME})
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${MAIN_OPENMM_LIB})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${SHARED_MBPOL_TARGET} optimized           ${SHARED_MBPOL_TARGET})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY")

INSTALL(TARGETS ${SHARED_TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/plugins)
SUBDIRS (tests)
