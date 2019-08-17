set(CMAKE_CXX_STANDARD 17) # -std=c++17

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall")

add_executable("simple-cooley-tukey01" ./simple-cooley-tukey01.cc)
#add_executable()
