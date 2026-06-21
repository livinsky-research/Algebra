file(GLOB EXAMPLE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

set(EXAMPLE_TARGET_NAMES
    Trees.cpp gen_trees
    RamseyExtension.cpp one_vertex
    RamseyGluing.cpp gluing
    Turan.cpp turan
)

foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_FILE ${EXAMPLE_SOURCE} NAME)
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)

    set(TARGET_NAME ${EXAMPLE_NAME})
    list(FIND EXAMPLE_TARGET_NAMES ${EXAMPLE_FILE} NAME_INDEX)
    if(NOT NAME_INDEX EQUAL -1)
        math(EXPR TARGET_INDEX "${NAME_INDEX} + 1")
        list(GET EXAMPLE_TARGET_NAMES ${TARGET_INDEX} TARGET_NAME)
    endif()

    add_executable(${TARGET_NAME} ${EXAMPLE_SOURCE})
    target_link_libraries(${TARGET_NAME} PRIVATE source)
    target_include_directories(${TARGET_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/inc)
endforeach()
