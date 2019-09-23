cmake_minimum_required(VERSION 3.2)
project(RSquare)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/main.cpp
        src/FileTypeCheck.h
        src/SummaryData.h src/SummaryData.cpp
        src/helperFunctions.h)
find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)

add_executable(RSquare ${SOURCE_FILES})

target_link_libraries(RSquare ${STATGEN_LIBRARY} ${ZLIB_LIBRARIES})

install(TARGETS RSquare RUNTIME DESTINATION bin)
