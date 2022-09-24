cmake_minimum_required(VERSION 3.8.2)

project(vectron)

set(SOURCE_FILES
    src/main.cpp
    src/renderer.cpp
    src/beam.cpp
    src/game.cpp)

set(CMAKE_CXX_STANDARD 11)

set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(SDL2 REQUIRED)

add_executable(main ${SOURCE_FILES})

if(APPLE)
    include_directories(${CMAKE_SOURCE_DIR} ${SDL2_INCLUDE_DIRS})
    target_link_libraries(main ${SDL2_LIBRARIES} OpenGL::GL GLEW::GLEW)
    set_target_properties(main PROPERTIES INSTALL_RPATH "@executable_path/")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH 1)
    #install(TARGETS main DESTINATION ${CMAKE_SOURCE_DIR})
else()
    target_link_libraries(main SDL2::SDL2 OpenGL::GL GLEW::GLEW)
endif()
