include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(simulator main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(simulator mrcimg ${OpenCV_LIBS})
