cmake_minimum_required(VERSION 3.10)
project(MultiphaseFlow)

set(CMAKE_CXX_STANDARD 11)

# 指定 C++ 标准

# 添加您的源代码文件
add_executable(MultiphaseFlow Multiphase.cpp common.h)


# 指定头文件目录
include_directories(include)

# 如果你的commen.h文件在include目录下，可以像下面这样指定头文件目录：
#include_directories(include)
# 查找线程库
find_package(Threads REQUIRED)

# 如果找到了线程库，则链接到你的可执行文件
if(THREADS_FOUND)
    target_link_libraries(MultiphaseFlow Threads::Threads)
endif()
