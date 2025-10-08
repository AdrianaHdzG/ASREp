# Define minimum cmake version
cmake_minimum_required( VERSION 3.12 )

# Define project name
project(ASRElib VERSION 1.0)

# add_definitions(-DPRINT_INT_RESULTS=1)

# Declare an library
add_library(ASRElibTimoBeam SHARED ASRElibTimoBeam.cpp)
add_library(ASRElib2DFrame SHARED ASRElib2DFrame.cpp)
add_library(ASRElib3DSolid SHARED ASRElib3DSolid.cpp)

# Include the headers in inc
target_include_directories(ASRElibTimoBeam PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")
target_include_directories(ASRElib2DFrame PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")
target_include_directories(ASRElib3DSolid PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")

# Create an INTERFACE target for the header-only library Eigen
add_library(Eigen INTERFACE)

# Specify the include directories for the header files of Eigen
target_include_directories(Eigen INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/Eigen")

# Link Eigen library to ASRElib
target_link_libraries(ASRElibTimoBeam PRIVATE Eigen)
target_link_libraries(ASRElib2DFrame PRIVATE Eigen)
target_link_libraries(ASRElib3DSolid PRIVATE Eigen)


get_filename_component(parent_dir "${CMAKE_CURRENT_SOURCE_DIR}" DIRECTORY)

# If Windows add a MACRO _MSVC
IF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")
    message("   ************************************")
    message("   Compiling on x64 Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _MSVC)
    target_compile_definitions(ASRElib2DFrame PRIVATE _MSVC)
    target_compile_definitions(ASRElib3DSolid PRIVATE _MSVC)
    # Define the dir to install the compiled dll
    install( TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/win32")
    install( TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/win32")
    install( TARGETS ASRElib3DSolid DESTINATION "${parent_dir}/bin/win32")
ELSEIF(MSVC)
    message("   ************************************")
    message("   Compiling on ${CMAKE_VS_PLATFORM_NAME} Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _MSVC)
    install( TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/win32")
    target_compile_definitions(ASRElib2DFrame PRIVATE _MSVC)
    install( TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/win32")
    target_compile_definitions(ASRElib3DSolid PRIVATE _MSVC)
    install( TARGETS ASRElib3DSolid DESTINATION "${parent_dir}/bin/win32")
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("   ************************************")
    message("   Compiling using GNU")
    message("   ************************************")
    target_compile_definitions(ASRElibTimoBeam PRIVATE _GCC)
    target_compile_definitions(ASRElib2DFrame PRIVATE _GCC)
    target_compile_definitions(ASRElib3DSolid PRIVATE _GCC)
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message(STATUS "Running on macOS")
        # if(CMAKE_C_COMPILER_TARGET MATCHES "x86_64")
        #     message(STATUS "Running on macOS x86_64")
        #     install( TARGETS ASRElib DESTINATION "${parent_dir}/bin/macOS")
        # ELSEIF(CMAKE_C_COMPILER_TARGET MATCHES "arm" OR CMAKE_C_COMPILER_TARGET MATCHES "aarch64")
        #     message(STATUS "Compiling for ARM target.")
        #     install( TARGETS ASRElib DESTINATION "${parent_dir}/bin/macOS_arm")
        # ELSE()
        #     message(STATUS "Unknown target architecture: ${CMAKE_C_COMPILER_TARGET}")
        # ENDIF()
        install( TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/macOS")
        install( TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/macOS")
        install( TARGETS ASRElib3DSolid DESTINATION "${parent_dir}/bin/macOS")
    ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        message(STATUS "Running on Linux")
        install( TARGETS ASRElibTimoBeam DESTINATION "${parent_dir}/bin/unix")
        install( TARGETS ASRElib2DFrame DESTINATION "${parent_dir}/bin/unix")
        install( TARGETS ASRElib3DSolid DESTINATION "${parent_dir}/bin/unix")
    ENDIF()
ELSE(TRUE)
    message("   ************************************")
    message("   Compiling on non-MSVC nor GNU")
    message("   ************************************")
ENDIF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")



