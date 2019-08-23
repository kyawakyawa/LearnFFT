set(CMAKE_CXX_STANDARD 17) # -std=c++17

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -Wall")

if (WITH_TCMALLOC)
  set(COMMON_LIBS tcmalloc)
endif ()

add_executable("bitrev01" ./bitrev01.cc)
add_executable("bitrev02" ./bitrev02.cc)
add_executable("number-functions01" ./number-functions01.cc)
add_executable("simple-cooley-tukey01" ./simple-cooley-tukey01.cc)
add_executable("simple-cooley-tukey02" ./simple-cooley-tukey02.cc)
add_executable("split-radix-index" ./split-radix-index.cc)
add_executable("transpose01" ./transpose01.cc)
add_executable("transpose02" ./transpose02.cc)

target_link_libraries("bitrev01"              ${COMMON_LIBS}) # Link Libraries
target_link_libraries("bitrev02"              ${COMMON_LIBS}) # Link Libraries
target_link_libraries("number-functions01"    ${COMMON_LIBS}) # Link Libraries
target_link_libraries("simple-cooley-tukey01" ${COMMON_LIBS}) # Link Libraries
target_link_libraries("simple-cooley-tukey02" ${COMMON_LIBS}) # Link Libraries
target_link_libraries("split-radix-index"     ${COMMON_LIBS}) # Link Libraries
target_link_libraries("transpose01"           ${COMMON_LIBS}) # Link Libraries
target_link_libraries("transpose02"           ${COMMON_LIBS}) # Link Libraries
#add_executable()
