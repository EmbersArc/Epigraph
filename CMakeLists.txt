cmake_minimum_required(VERSION 3.5)
project(epigraph)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED NO_MODULE)

set(TARGET_INCLUDE include solvers/wrappers/include)

set(SOCP_SOURCES
    src/parameter.cpp
    src/variable.cpp
    src/expressions.cpp
    src/constraint.cpp
    src/problem.cpp
    solvers/wrappers/src/socpWrapperBase.cpp
    solvers/wrappers/src/qpWrapperBase.cpp
    solvers/wrappers/src/osqpWrapper.cpp
    solvers/wrappers/src/ecosWrapper.cpp
    solvers/wrappers/src/eicosWrapper.cpp)

# ==== Solvers ====
# Make sure the repository is cloned with --recurse-submodules.

# Note that enabling certain solvers automatically alters the license.
option(ENABLE_OSQP "Enable OSQP solver support. Apache-2.0 License." FALSE)
option(ENABLE_ECOS "Enable ECOS solver support. GPLv3 License." FALSE)
option(ENABLE_EiCOS "Enable EiCOS solver support. MIT License." FALSE)

add_library(epigraph SHARED ${SOCP_SOURCES})

set(TARGET_INCLUDE_SYS solvers)

# OSQP
if(ENABLE_OSQP)
  add_subdirectory(solvers/osqp)
  target_link_libraries(epigraph osqp)
  add_compile_definitions(ENABLE_OSQP)
  message("Epigraph: OSQP enabled.")
endif()

# ECOS
if(ENABLE_ECOS)
  execute_process(COMMAND make shared
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/solvers/ecos)
  list(APPEND TARGET_INCLUDE_SYS solvers/ecos/external/SuiteSparse_config)
  target_link_libraries(epigraph
                        ${CMAKE_CURRENT_LIST_DIR}/solvers/ecos/libecos.so)
  add_compile_definitions(ENABLE_ECOS)
  message("Epigraph: ECOS enabled.")
endif()

# EiCOS
if(ENABLE_EICOS)
  add_subdirectory(solvers/eicos)
  target_link_libraries(epigraph eicos)
  add_compile_definitions(ENABLE_EICOS)
  message("Epigraph: EiCOS enabled.")
endif()

target_include_directories(epigraph PUBLIC ${TARGET_INCLUDE})
target_include_directories(epigraph SYSTEM PUBLIC ${TARGET_INCLUDE_SYS})

set(DEBUG_OPTIONS -Wall -Wextra -Wpedantic)
set(RELEASE_OPTIONS -O2)

target_link_libraries(epigraph Eigen3::Eigen)

# ==== Testing ====
find_package(Catch2)
option(ENABLE_TESTING "Build tests. Requires Catch2 and fmt to be installed."
       FALSE)
if(Catch2_FOUND AND ENABLE_TESTING)
  if(NOT
     (ENABLE_OSQP
      AND ENABLE_ECOS
      AND ENABLE_EICOS))
    message(FATAL_ERROR "Epigraph: Enable all solvers to enable testing.")
  endif()

  include_directories(epigraph tests)
  add_executable(run_tests tests/run_tests.cpp)
  target_link_libraries(run_tests epigraph)
  message("Epigraph: Testing is enabled.")

  find_package(fmt REQUIRED)
  if(${fmt_FOUND})
    message(STATUS "Epigraph: Found fmt.")
    target_link_libraries(run_tests fmt::fmt)
  endif()

  # Coverage
  if(CMAKE_COMPILER_IS_GNUCXX)
    option(ENABLE_COVERAGE "Enable coverage reporting for gcc/clang" FALSE)
  endif()

  if(ENABLE_COVERAGE)
    message("Epigraph: Code coverage is enabled.")
    list(APPEND DEBUG_OPTIONS --coverage -O0)
  endif()

  target_link_libraries(epigraph --coverage)
endif()
# ====

target_compile_options(epigraph PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(epigraph
                       PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
