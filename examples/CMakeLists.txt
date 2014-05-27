CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# set OpenMM install folder
SET(OpenMM_INSTALL_DIR "/usr/local/openmm")

INCLUDE_DIRECTORIES (${OpenMM_INSTALL_DIR}/include/) 

ADD_LIBRARY( openmmlib SHARED IMPORTED )
SET_TARGET_PROPERTIES( openmmlib PROPERTIES IMPORTED_LOCATION (${OpenMM_INSTALL_DIR}/lib/ ) )

ADD_EXECUTABLE (Simulate14WaterCluster Simulate14WaterCluster.cpp)

TARGET_LINK_LIBRARIES(Simulate14WaterCluster openmmlib)
