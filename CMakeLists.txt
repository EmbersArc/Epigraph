cmake_minimum_required(VERSION 3.5)
project(epigraph)

set(CMAKE_CXX_STANDARD 17)

add_library(epigraph SHARED)
target_include_directories(epigraph PUBLIC include)
target_sources(epigraph PRIVATE
    src/parameter.cpp
    src/variable.cpp
    src/expressions.cpp
    src/constraint.cpp
    src/problem.cpp

    src/wrappers/wrapperBase.cpp
    src/wrappers/socpWrapperBase.cpp
    src/wrappers/qpWrapperBase.cpp
)

# ==== Solvers ====
# Make sure the repository is cloned with --recurse-submodules.

# Note that enabling certain solvers automatically alters the license.
# Set ENABLE_OSQP to TRUE to enable OSQP solver support. Apache-2.0 License.
# Set ENABLE_ECOS to TRUE to enable ECOS solver support. GPLv3 License.
find_package(Git QUIET)

if(NOT ENABLE_OSQP AND NOT ENABLE_ECOS)
  message(FATAL_ERROR "Epigraph: No solvers enabled. Set the appropriate variables for the desired solvers and re-run cmake.")
endif()

# OSQP
if(ENABLE_OSQP)
  find_package(osqp QUIET)
  if(osqp_FOUND)
    message("Using installed OSQP.")
    target_link_libraries(epigraph osqp::osqp)
  else()
    message("Using OSQP from subdirectory.")
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/solvers/osqp
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to load the OSQP submodule with ${GIT_SUBMOD_RESULT}.")
      endif()
      add_subdirectory(solvers/osqp)
      target_link_libraries(epigraph osqp)
    else()
      message(FATAL_ERROR "Git is not available or configured incorrectly.")
    endif()
  endif()

  target_sources(epigraph PRIVATE src/wrappers/osqpWrapper.cpp)
  target_compile_definitions(epigraph PUBLIC ENABLE_OSQP)
  message("Epigraph: OSQP enabled.")
else()
  message("Epigraph: OSQP disabled.")
endif()

# ECOS
if(ENABLE_ECOS)
  find_package(ecos QUIET)
  if(ecos_FOUND)
    message("Using installed ECOS.")
    target_link_libraries(epigraph ecos::ecos)
  else()
    message("Using ECOS from subdirectory.")
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/solvers/ecos
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to load the ECOS submodule with ${GIT_SUBMOD_RESULT}.")
      endif()
      add_subdirectory(solvers/ecos)
      target_link_libraries(epigraph ecos)
    else()
      message(FATAL_ERROR "Git is not available or configured incorrectly.")
    endif()
  endif()

  target_sources(epigraph PRIVATE src/wrappers/ecosWrapper.cpp)
  target_compile_definitions(epigraph PUBLIC ENABLE_ECOS)
  message("Epigraph: ECOS enabled.")
else()
  message("Epigraph: ECOS disabled.")
endif()

# Set compiler flags
set(DEBUG_OPTIONS "")
set(RELEASE_OPTIONS "")

if(MSVC)
    list(APPEND DEBUG_OPTIONS "/Od" "/Wall")
    list(APPEND RELEASE_OPTIONS "/O2")
else()
    list(APPEND DEBUG_OPTIONS "-O0" "-Wall" "-Wextra")
    list(APPEND RELEASE_OPTIONS "-O3")
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(epigraph Eigen3::Eigen)

# ==== Testing ====
find_package(Catch2)
option(ENABLE_TESTING "Build tests. Requires Catch2 and fmt to be installed."
       FALSE)
if(Catch2_FOUND AND ENABLE_TESTING)
  if(NOT (ENABLE_OSQP AND ENABLE_ECOS))
    message(FATAL_ERROR "Epigraph: Enable all solvers to enable testing.")
  endif()

  include_directories(epigraph tests)
  add_executable(run_tests tests/run_tests.cpp)
  target_link_libraries(run_tests epigraph)
  message("Epigraph: Testing is enabled.")

  # fmt
  find_package(fmt REQUIRED)
  if(${fmt_FOUND})
    message(STATUS "Epigraph: Found fmt.")
    target_link_libraries(run_tests fmt::fmt)
  endif()

  # Coverage
  if(CMAKE_COMPILER_IS_GNUCXX)
    option(ENABLE_COVERAGE "Enable coverage reporting for gcc/clang" FALSE)
    if(ENABLE_COVERAGE)
      message("Epigraph: Code coverage is enabled.")
      list(APPEND DEBUG_OPTIONS --coverage -O0)
      target_link_libraries(epigraph --coverage)
    endif()
  endif()
endif()
# ====

target_compile_options(epigraph PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(epigraph PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
