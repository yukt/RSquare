cmake_minimum_required(VERSION 3.2)
project(Rsquare)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/main.cpp
        src/RSquare.h src/RSquare.cpp
        src/Haplotype.h src/Haplotype.cpp)
find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)

add_executable(Rsquare ${SOURCE_FILES})

target_link_libraries(Rsquare ${STATGEN_LIBRARY} ${ZLIB_LIBRARIES})

install(TARGETS Rsquare RUNTIME DESTINATION bin)
