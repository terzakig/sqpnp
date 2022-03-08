#project(RobustExample)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++14)

#set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}) # RansacLib

set(sources main.cpp robust_pose_pnp.cpp )

# Append list of libs with SQPnP (static library)
list( APPEND LIBS SQPnP )

add_executable(robust_sqpnp_example ${sources})

target_link_libraries(robust_sqpnp_example ${LIBS} )
