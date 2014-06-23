cmake_minimum_required(VERSION 2.8.3)
project(kr_attitude_eskf)
find_package(catkin REQUIRED COMPONENTS roscpp tf geometry_msgs sensor_msgs
                                        message_filters message_generation
                                        eigen_conversions kr_math)

# generate messages
add_service_files(FILES BeginCalibration.srv)
add_message_files(FILES Status.msg)
generate_messages(DEPENDENCIES std_msgs)

catkin_package()

# include boost
find_package(Boost REQUIRED)
link_directories(${Boost_LIBRARY_DIR})
include_directories(${Boost_INCLUDE_DIR})

# include headers
include_directories(${kr_math_INCLUDE_DIRS})

# build shared library
add_library(krAttitudeESKF src/AttitudeESKF.cpp)

# build main filter
add_executable(${PROJECT_NAME} src/kr_attitude_eskf.cpp)

target_link_libraries(${PROJECT_NAME} krAttitudeESKF
                                      ${roscpp_LIBRARIES}
                                      ${galt_common_LIBRARIES}
                                      ${tf_LIBRARIES}
                                      ${eigen_conversions_LIBRARIES}
                                      ${message_filters_LIBRARIES})

# ensures messages are generated before hand
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS})
