cmake_minimum_required(VERSION 3.0.0)
project(PPMViewer VERSION 0.1.0)

find_package(OpenGL REQUIRED)

set(LIBS
    ${OPENGL_LIBRARY}
    glfw
    ${GLFW_LIBRARIES}
)

include_directories(external/glad/include)

add_subdirectory(external)

add_executable(PPMViewer ppmviewer.cpp)
target_link_libraries(PPMViewer ${LIBS})
