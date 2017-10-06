cmake_minimum_required(VERSION 3.2)
project(Rsquare)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/main.cpp
        src/FileTypeCheck.h
        src/SummaryData.h src/SummaryData.cpp
        src/helperFunctions.h)
find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)

add_executable(Rsquare ${SOURCE_FILES})

target_link_libraries(Rsquare ${STATGEN_LIBRARY} ${ZLIB_LIBRARIES})

install(TARGETS Rsquare RUNTIME DESTINATION bin)
