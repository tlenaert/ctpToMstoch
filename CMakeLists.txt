cmake_minimum_required(VERSION 3.16.0)
project (ctpstoch)

list(APPEND flags "-march=native")

set(SOURCE_FILES main.cpp
		 kernel.cpp
                 strategy.cpp
                 ctpgame.cpp
                 ctpdata.cpp
                 rangen.cpp)

add_executable(ctpstoch ${SOURCE_FILES})

target_compile_options(ctpstoch
        PRIVATE
                ${flags}
)

find_package(GSL REQUIRED)    
target_link_libraries(ctpstoch GSL::gsl GSL::gslcblas)
