include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(mask main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(mask mrcimg ${OpenCV_LIBS})
