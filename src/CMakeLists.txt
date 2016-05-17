########## add subdirectories ############

set(APP_SRC 
    ${APP_SRC}
    model.cpp    
    main.cpp
    Splitter.cpp
    spatial_hash.cpp
    unfolder.cpp
    OverlappingChecker.cpp
)

set(APP_NAME unfolder)

add_subdirectory(GL)
add_subdirectory(mathtool)
add_subdirectory(modelgraph)
add_subdirectory(util)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

include_directories(    
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/src/GL"
    "${CMAKE_SOURCE_DIR}/src/mathtool"
    "${CMAKE_SOURCE_DIR}/src/modelgraph" 
    "${CMAKE_SOURCE_DIR}/libs")

link_directories(    
    "${CMAKE_SOURCE_DIR}/libs"
    "${CMAKE_SOURCE_DIR}/libs/GL"
    "${CMAKE_SOURCE_DIR}/libs/itree"
)

add_executable(${APP_NAME} ${APP_SRC})

target_link_libraries(${APP_NAME}
                    itree
                    ${OPENGL_LIBRARIES}
                    ${GLUT_LIBRARY})

if(WIN32)
    target_link_libraries(${APP_NAME} glut32.lib)
endif()
