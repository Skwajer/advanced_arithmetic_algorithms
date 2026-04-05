cmake_minimum_required(VERSION 3.16)
project(monomial_ideal)

set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

qt_add_executable(monomial_ideal
    main.cpp
    mainwindow.cpp
    mainwindow.hpp
)

target_link_libraries(monomial_ideal PRIVATE Qt6::Widgets)
