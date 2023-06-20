cmake_minimum_required(VERSION 3.10)
project(AS4PM CXX)

set(SOURCES
        main.cpp
        Particle.cpp
        matplotlib-cpp/matplotlibcpp.h
        )
        #test.cpp

add_executable(AS4PM ${SOURCES})

target_compile_features(AS4PM PRIVATE cxx_std_11)

# Add the include directories for Python 3.10
target_include_directories(AS4PM PRIVATE
        /opt/homebrew/Cellar/libomp/16.0.2/include
        /opt/homebrew/Cellar/python@3.9/3.9.16/Frameworks/Python.framework/Versions/3.9/include/python3.9
        /opt/homebrew/lib/python3.9/site-packages/numpy/core/include
        )

# Add the library directories for Python 3.10
target_link_directories(AS4PM PRIVATE
        /opt/homebrew/Cellar/python@3.9/3.9.16/Frameworks/Python.framework/Versions/3.9/lib
        )

# Link against the Python 3.10 library
target_link_libraries(AS4PM PRIVATE
        python3.9

        )

set_target_properties(AS4PM PROPERTIES OUTPUT_NAME "AS4PM")

target_compile_options(AS4PM PRIVATE -Wall -Wextra -pedantic)
