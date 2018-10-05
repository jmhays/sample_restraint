# Defines targets for the C++ restraints implemented here. These CMake targets are used by the
# unit tests and by the Python module target defined in ../pythonmodule/CMakeLists.txt

# Create a shared object library for a simple restraint.
add_library(harmonicpotential STATIC
            harmonicpotential.h
            harmonicpotential.cpp
            sessionresources.cpp)
set_target_properties(harmonicpotential PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(harmonicpotential PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(harmonicpotential PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(harmonicpotential PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(harmonicpotential PRIVATE Gromacs::gmxapi)


# Create a shared object library for our restrained ensemble plugin.
add_library(ensemblepotential STATIC
            ensemblepotential.h
            ensemblepotential.cpp
            sessionresources.cpp)
set_target_properties(ensemblepotential PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(ensemblepotential PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(ensemblepotential PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(ensemblepotential PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(ensemblepotential PRIVATE Gromacs::gmxapi)


# Add metadynamics code
add_library(ebmetad STATIC
            ebmetad.h
            ebmetad.cpp
            sessionresources.cpp)
set_target_properties(ebmetad PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(ebmetad PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(ebmetad PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(ebmetad PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(ebmetad PRIVATE Gromacs::gmxapi)

# Add linear potential
add_library(linear STATIC
            linearpotential.h
            linearpotential.cpp
            sessionresources.cpp)
set_target_properties(linear PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(linear PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(linear PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(linear PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(linear PRIVATE Gromacs::gmxapi)

# Add linearstop potential
add_library(linearstop STATIC
            linearstoppotential.h
            linearstoppotential.cpp
            sessionresources.cpp)
set_target_properties(linearstop PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(linearstop PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(linearstop PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(linearstop PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(linearstop PRIVATE Gromacs::gmxapi)


# Add BRER code
add_library(brer STATIC
            brerpotential.h
            brerpotential.cpp
            sessionresources.cpp)
set_target_properties(brer PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(brer PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(brer PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(brer PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(brer PRIVATE Gromacs::gmxapi)
