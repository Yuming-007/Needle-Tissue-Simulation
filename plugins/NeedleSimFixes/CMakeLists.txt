cmake_minimum_required(VERSION 3.22)
project(NeedleSimFixes VERSION 0.1 LANGUAGES CXX)

# 针对 SOFA v25.12 已知 bug 的最小修补。不修改 SOFA 源码，只通过继承覆盖有问题的函数。
find_package(Sofa.Component.SolidMechanics.FEM.Elastic REQUIRED)

add_library(${PROJECT_NAME} SHARED
    src/NeedleSimFixes/init.cpp
    src/NeedleSimFixes/FastTetrahedralCorotationalForceFieldFixed.cpp)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
target_include_directories(${PROJECT_NAME} PUBLIC src)
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Component.SolidMechanics.FEM.Elastic)
