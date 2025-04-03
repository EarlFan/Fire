cmake_minimum_required(VERSION 3.10)
project(Fire)

# define the total species number of the gasous mixture in the simulation
set(TOTAL_SPECIES_NUMBER 10)

# --------------------------------------------
# Automatic Cantera Environment Detection
# --------------------------------------------

# 1. Mandatory conda environment query configuration
set(CANTERA_ENV_NAME "ct" CACHE STRING "Name of conda environment containing Cantera")

# 2. Locate conda executable
find_program(CONDA_EXE conda REQUIRED)
message(STATUS "Found conda: ${CONDA_EXE}")

# 3. Query conda environments
message(STATUS "Searching for conda env: ${CANTERA_ENV_NAME}")
execute_process(
  COMMAND bash -c "${CONDA_EXE} info --envs | grep -E '^ct(\\s+|\\*\\s+)' | awk -F'[*[:space:]]+' '/ct/ {print $(NF)}'"
  OUTPUT_VARIABLE CANTERA_ENV_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE retcode
)

# 4. Error handling for environment query
if(NOT retcode EQUAL 0)
    message(FATAL_ERROR "Failed to find conda environment '${CANTERA_ENV_NAME}'. Debug steps:
    1. Run 'conda info --envs' manually to verify environment existence
    2. Check conda executable path: which conda")
endif()

if("${CANTERA_ENV_PATH}" STREQUAL "")
    message(FATAL_ERROR "Conda environment '${CANTERA_ENV_NAME}' path is empty!
    Possible causes:
    - Environment not created through 'conda create'
    - Environment metadata corruption")
endif()

# 5. Validate path integrity
if(NOT IS_DIRECTORY "${CANTERA_ENV_PATH}")
    message(FATAL_ERROR "Invalid conda environment path: ${CANTERA_ENV_PATH}
    Verify environment installation status")
endif()
message(STATUS "Cantera environment path: ${CANTERA_ENV_PATH}")

add_compile_definitions(NS=${TOTAL_SPECIES_NUMBER})

# Generate config file for shell scripts
file(WRITE ${CMAKE_SOURCE_DIR}/config.sh "export NS=${TOTAL_SPECIES_NUMBER}\n")


# Set the Cantera include directory
set(CANTERA_INCLUDE ${CANTERA_ENV_PATH}/include)
set(CANTERA_LIB ${CANTERA_ENV_PATH}/lib)

# Set the source file extension and paths
set(SRC_EXT cpp)
set(SRC_PATH ../../src)
set(BUILD_PATH build)
set(BIN_PATH ${BUILD_PATH}/bin)

# Set the build type (default to Release if not specified)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set the compiler and compiler flags
set(CMAKE_CXX_COMPILER mpic++)

# Set compiler flags based on the build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --debug -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe -isystem ${CANTERA_INCLUDE} -pthread -std=c++17 -Wno-c++17-extensions -Wno-return-type")
  set(EXE_NAME Fire_Debug)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe -isystem ${CANTERA_INCLUDE} -pthread -std=c++17 -Wno-c++17-extensions -Wno-return-type")
  set(EXE_NAME Fire_Release)
endif()

# Set the executable name
set(BIN_NAME "${EXE_NAME}_${TOTAL_SPECIES_NUMBER}Species")

# Include cantera directories
include_directories(${CANTERA_INCLUDE})
link_directories(${CANTERA_LIB})

# Find all source files in the source directory
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "${SRC_PATH}/*.${SRC_EXT}")

# Set the output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${BIN_PATH})

# Add the executable
add_executable(${BIN_NAME} ${SOURCES})

# Set the target properties
set_target_properties(${BIN_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${BIN_PATH}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/${BIN_PATH}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/${BIN_PATH}
)

# Add custom commands for clean and copy
add_custom_target(clean-target
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/${BUILD_PATH}
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${BIN_NAME}
)

add_custom_command(TARGET ${BIN_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/${BIN_PATH}/${BIN_NAME} ${CMAKE_BINARY_DIR}/${BIN_NAME}
)

# Create necessary directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${BUILD_PATH})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${BIN_PATH})

target_link_libraries(${BIN_NAME} cantera dl)