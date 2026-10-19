cmake_minimum_required(VERSION 3.20)
project(browser_qchem VERSION 0.2.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/libcint/CMakeLists.txt")
  message(FATAL_ERROR "Run: git submodule update --init --recursive")
endif()

# Keep the upstream sources unmodified. Use portable scalar C for the first port.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_MARCH_NATIVE OFF CACHE BOOL "" FORCE)
set(WITH_FORTRAN OFF CACHE BOOL "" FORCE)
set(WITH_CINT2_INTERFACE OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLE OFF CACHE BOOL "" FORCE)
set(ENABLE_TEST OFF CACHE BOOL "" FORCE)
if(EMSCRIPTEN)
  set(CMAKE_DISABLE_FIND_PACKAGE_QUADMATH TRUE)
endif()
add_subdirectory(vendor/libcint EXCLUDE_FROM_ALL)
set_target_properties(cint PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(EMSCRIPTEN)
  add_executable(qc_core src/overlap.c src/integrals.c)
  set_target_properties(qc_core PROPERTIES
    SUFFIX ".mjs"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/dist")
  target_link_options(qc_core PRIVATE
    --no-entry
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sEXPORT_NAME=createQChem"
    "-sENVIRONMENT=web,worker,node"
    "-sEXPORTED_FUNCTIONS=['_h2_overlap','_qc_compute','_qc_error','_qc_nao','_qc_enuc','_qc_data','_qc_free','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['HEAP32','HEAPF64','UTF8ToString']"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sSTACK_SIZE=1048576")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/dist")
  foreach(asset index.html app.mjs worker.mjs molecule.mjs basis.mjs engine.mjs)
    configure_file("web/${asset}" "${CMAKE_BINARY_DIR}/dist/${asset}" COPYONLY)
  endforeach()
else()
  add_library(qc_core SHARED src/overlap.c src/integrals.c)
endif()

target_include_directories(qc_core PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/vendor/libcint/include"
  "${CMAKE_CURRENT_BINARY_DIR}/vendor/libcint/include")
target_link_libraries(qc_core PRIVATE cint m)
