add_library(SQPnP sqpnp.cpp)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

#list( APPEND CMAKE_MODULE_PATH
#  ${CMAKE_SOURCE_DIR}/cmake
#)

find_package(Eigen3 REQUIRED)
list(APPEND INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

include_directories(${INCLUDE_DIRS})

# when other libraries or executables link to SQPnP
target_include_directories(SQPnP PUBLIC ${INCLUDE_DIRS})
