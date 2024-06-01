cmake_minimum_required(VERSION 3.10)
project(Fire)

# define the path to cantera envs of conda
set(CANTERA_ENV_PATH /home/fane/anaconda3/envs/ct-test)
# define the total species number of the gasous mixture in the simulation
set(TOTAL_SPECIES_NUMBER 9)

add_compile_definitions(NS=${TOTAL_SPECIES_NUMBER})

# Set the Cantera include directory
set(CANTERA_INCLUDE ${CANTERA_ENV_PATH}/include)
set(CANTERA_LIB ${CANTERA_ENV_PATH}/lib)

# Set the source file extension and paths
set(SRC_EXT cpp)
set(SRC_PATH src)
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

# include_directories("/home/fane/anaconda3/envs/ct-cf-dev/include")

target_link_libraries(${BIN_NAME} cantera dl)