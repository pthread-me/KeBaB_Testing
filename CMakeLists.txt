cmake_minimum_required(VERSION 4.0)

project(AES)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE Debug)

if(False)
set(SRCs src/nt_hash.cpp)
add_library(AES STATIC ${SRCs})
target_include_directories(AES PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(AES PUBLIC 
	-Wall
	-Wextra
	-Wconversion
	-Waddress
	-Wpedantic
	-O0

	# required
	-march=native
	-maes
	-msse4
	-mvaes

	#asm
)
endif()

# for testing cause why would i ever use ctest :|
add_executable(main src/main.cpp)
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(main PUBLIC AES)
target_compile_options(main PUBLIC
	-save-temps
	-fverbose-asm
	-march=native
	-msse4
  -O0
)

