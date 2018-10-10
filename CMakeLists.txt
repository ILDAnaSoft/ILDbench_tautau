########################################################
# cmake file for building ILDbench_tautau analysis package
# @author Ete Remi, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
PROJECT( ILDbench_tautau )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################
FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( DD4hep )
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} )
ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )

FIND_PACKAGE( MarlinUtil 1.0 REQUIRED ) # minimum required MarlinUtil version
INCLUDE_DIRECTORIES( SYSTEM ${MarlinUtil_INCLUDE_DIRS} )
LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )
ADD_DEFINITIONS( ${MarlinUtil_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED ) 
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

FIND_PACKAGE( LCFIVertex REQUIRED )
INCLUDE_DIRECTORIES( ${LCFIVertex_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${LCFIVertex_ROOT}/vertex_lcfi ${LCFIVertex_ROOT}/boost )
LINK_LIBRARIES( ${LCFIVertex_LIBRARIES} )

### tauola stuff
INCLUDE_DIRECTORIES( /home/ilc/jeans/LHAPDF/lhapdf-5.9.1/include )
INCLUDE_DIRECTORIES( /home/ilc/jeans/tauola/v116c/TAUOLA/include )
LINK_DIRECTORIES( /home/ilc/jeans/tauola/v116c/TAUOLA/lib/  /home/ilc/jeans/LHAPDF/lhapdf-5.9.1/install/lib )
LINK_LIBRARIES( LHAPDF TauolaCxxInterface TauolaFortran TauolaTauSpinner )



# Compile you project
ADD_SUBDIRECTORY( source )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
