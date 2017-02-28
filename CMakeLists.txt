CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project("math" C)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -std=c99 -Wall -Werror")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c99 -O2 -ffast-math -funroll-loops")

if(WIN32)
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif(WIN32)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message("Default build type is RELEASE!")
endif(NOT CMAKE_BUILD_TYPE)

file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.c)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME} ${SRC})
