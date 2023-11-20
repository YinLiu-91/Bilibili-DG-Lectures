# CMake 版本
cmake_minimum_required(VERSION 3.16.3)
# 指定Fortran编译器
# set(CMAKE_Fortran_COMPILER "/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
# 指定语言
enable_language(Fortran)
# # MKL库的位置
set(mkllib "/opt/intel/oneapi/mkl/2022.0.2/lib/intel64")
# #链接该库位置
link_directories(${mkllib})
#不同的编译器设置
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # set(dialect "-ffree-form -std=f2003 -fimplicit-none")
    set(dialect "-ffree-form -std=f2008 -ffree-line-length-none -freal-4-real-8 -mcmodel=medium")
    # set(bounds "-fbounds-check")
endif()
#本例所用ifort编译器设置 使用默认设置
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message("HERE")
    set(dialect "-stand f08 -free -implicitnone") 
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# 项目名
project(Example)
# 编译为debug模式
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

# 添加工作目录下子目录src中的源文件到SRC_FILES变量
file(GLOB  SRC_FILES
    "${PROJECT_SOURCE_DIR}/*.f90"
    # "${PROJECT_SOURCE_DIR}/src/*.for"
    )

# 生成可执行文件 
#${CMAKE_PROJECT_NAME}指项目名 Example
add_executable(${CMAKE_PROJECT_NAME}  ${SRC_FILES})
#链接到MKL库
# 注意没有前面的lib 因为编译器自动加 -l 表示lib
target_link_libraries(${CMAKE_PROJECT_NAME}                       
                        mkl_core
                        mkl_intel_lp64
                        mkl_sequential )  