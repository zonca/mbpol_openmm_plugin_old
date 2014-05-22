#set(GCCXML_EXTRA_ARGS "" CACHE STRING "Additional arguments to gccxml, such as '--gccxml-compiler;msvc8'")
#SET(GCCXML_ARGS) # start empty
#FOREACH(subdir ${API_INCLUDE_DIRS})
#    SET(GCCXML_ARGS ${GCCXML_ARGS} -I${subdir})
#ENDFOREACH(subdir)
#SET(GCCXML_ARGS ${GCCXML_ARGS} ${GCCXML_EXTRA_ARGS})

SET(SAXON_JAR ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar)

# OpenMMMBPol

ADD_CUSTOM_COMMAND(OUTPUT MBPolOpenMMApi.xml COMMAND ${GCCXML_PATH} ${GCCXML_ARGS} -I${CMAKE_SOURCE_DIR}/openmmapi/include -I${CMAKE_SOURCE_DIR}/olla/include ${CMAKE_SOURCE_DIR}/plugins/mbpol/openmmapi/include/OpenMMMBPol.h -fxml=MBPolOpenMMApi.xml)

ADD_CUSTOM_COMMAND(OUTPUT MBPolOpenMMCWrapper.h COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:MBPolOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/mbpol/wrappers/CWrapper_Header.xslt -o:MBPolOpenMMCWrapper.h DEPENDS MBPolOpenMMApi.xml)
#
ADD_CUSTOM_COMMAND(OUTPUT MBPolOpenMMCWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:MBPolOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/mbpol/wrappers/CWrapper_Source.xslt -o:MBPolOpenMMCWrapper.cpp DEPENDS MBPolOpenMMApi.xml)

ADD_CUSTOM_COMMAND(OUTPUT MBPolOpenMMFortranModule.f90 COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:MBPolOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/mbpol/wrappers/FortranWrapper_Header.xslt -o:MBPolOpenMMFortranModule.f90 DEPENDS MBPolOpenMMApi.xml ${CMAKE_SOURCE_DIR}/plugins/mbpol/wrappers/FortranWrapper_Header.xslt)

ADD_CUSTOM_COMMAND(OUTPUT MBPolOpenMMFortranWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:MBPolOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/mbpol/wrappers/FortranWrapper_Source.xslt -o:MBPolOpenMMFortranWrapper.cpp DEPENDS MBPolOpenMMApi.xml)

ADD_CUSTOM_TARGET(MBPolApiWrappers DEPENDS MBPolOpenMMCWrapper.h MBPolOpenMMCWrapper.cpp  MBPolOpenMMFortranModule.f90 MBPolOpenMMFortranWrapper.cpp)

INSTALL_FILES(/include FILES MBPolOpenMMCWrapper.h MBPolOpenMMFortranModule.f90 )
