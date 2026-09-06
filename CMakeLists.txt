cmake_minimum_required(VERSION 4.3)
project(Code)

set(CMAKE_CXX_STANDARD 14)

add_executable(Code
        .idea/.gitignore
        .idea/workspace.xml
        main.c
        shared_euler.c
        shared_euler.h
        eulerian_1d.h
        eulerian_1d.c
        array_handling.c
        array_handling.h
        io.c
        io.h)
