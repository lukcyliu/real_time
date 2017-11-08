cmake_minimum_required(VERSION 3.7)
project(real_time)

set(CMAKE_C_STANDARD 11)
set(PROJECT_NAME real_time)
set(SOURCE_FILES main.c NewKalmanFilter.c MahonyAHRS.c matrix_lc.c Queue.c Turnning.c)
add_executable(real_time ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} -lm -lpthread)