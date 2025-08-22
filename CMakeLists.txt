cmake_minimum_required(VERSION 3.20)

# ---- Project ----
project(AMSC-ExamProject
  VERSION 1.0
  LANGUAGES CXX
)

# ---- Options ----
option(USE_OMP        "Enable OpenMP (sets -DUSE_OMP if found)" ON)
option(ENABLE_NATIVE  "Enable -march=native (or equivalent) on non-MSVC compilers" OFF)
option(ENABLE_WARNINGS "Enable extra warnings" ON)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Sources / Includes ----
set(PROJECT_SOURCES
  src/mc_explicit.cpp
  src/main.cpp
)
add_executable(pricer ${PROJECT_SOURCES})
target_include_directories(pricer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ---- Warnings (portable) ----
if(ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(pricer PRIVATE /W4 /permissive-)
  else()
    target_compile_options(pricer PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wno-unused-parameter)
  endif()
endif()

# ---- march=native (optional) ----
if(ENABLE_NATIVE AND NOT MSVC)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
  if(HAS_MARCH_NATIVE)
    target_compile_options(pricer PRIVATE -march=native)
  endif()
endif()

# ---- OpenMP (portable) ----
if(USE_OMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(pricer PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(pricer PRIVATE USE_OMP)
  else()
    message(STATUS "OpenMP not found; building without OpenMP.")
  endif()
endif()

# ---- Release flags (keep CMake defaults but ensure -O3 for GCC/Clang) ----
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(pricer PRIVATE
    $<$<CONFIG:Release>:-O3>
  )
endif()

# ---- Clean all (custom target) ----
add_custom_target(clean_all
  COMMAND ${CMAKE_COMMAND} -E echo "Cleaning build artifacts and generated files..."
  # remove the binary if created in source root
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/pricer
  # remove only CSV/log/data outputs, not all *.txt
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/*.csv
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/*.dat
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/*.log
  # wipe the build directory
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/*
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# ---- Install target (optional) ----
include(GNUInstallDirs)
install(TARGETS pricer RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ---- Nice messages ----
message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "USE_OMP: ${USE_OMP}")
message(STATUS "ENABLE_NATIVE: ${ENABLE_NATIVE}")
message(STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")