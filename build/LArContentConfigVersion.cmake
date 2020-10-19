##############################################################################
# this file is parsed when FIND_PACKAGE is called with version argument
#
# @author Jan Engels, Desy IT
##############################################################################


SET( ${PACKAGE_FIND_NAME}_VERSION_MAJOR 03 )
SET( ${PACKAGE_FIND_NAME}_VERSION_MINOR 15 )
SET( ${PACKAGE_FIND_NAME}_VERSION_PATCH 08 )


INCLUDE( "/usera/afm67/2020/May/CosmicJorisMay2020/PandoraPFA/cmakemodules/MacroCheckPackageVersion.cmake" )
CHECK_PACKAGE_VERSION( ${PACKAGE_FIND_NAME} 03.15.08 )

