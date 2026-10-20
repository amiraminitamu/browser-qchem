cmake_minimum_required(VERSION 3.20)
project(browser_qchem VERSION 0.5.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/libcint/CMakeLists.txt")
  message(FATAL_ERROR "Run: git submodule update --init --recursive")
endif()
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen/Eigen/Core")
  message(FATAL_ERROR "Run: git submodule update --init --recursive (Eigen is required)")
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
  add_executable(qc_core src/overlap.c src/integrals.c src/rhf.cpp src/uhf.cpp src/mp2.cpp)
  set_target_properties(qc_core PROPERTIES
    SUFFIX ".mjs"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/dist")
  target_link_options(qc_core PRIVATE
    --no-entry
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sEXPORT_NAME=createQChem"
    "-sENVIRONMENT=web,worker,node"
    "-sEXPORTED_FUNCTIONS=['_h2_overlap','_qc_compute','_qc_error','_qc_nao','_qc_enuc','_qc_data','_qc_free','_qc_rhf_create','_qc_rhf_error','_qc_rhf_step','_qc_rhf_status','_qc_rhf_iteration','_qc_rhf_scalar','_qc_rhf_data','_qc_rhf_free','_qc_uhf_create','_qc_uhf_error','_qc_uhf_step','_qc_uhf_status','_qc_uhf_iteration','_qc_uhf_scalar','_qc_uhf_data','_qc_uhf_free','_qc_mp2_create','_qc_mp2_error','_qc_mp2_step','_qc_mp2_status','_qc_mp2_iteration','_qc_mp2_nocc','_qc_mp2_nvirt','_qc_mp2_scalar','_qc_mp2_free','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['HEAP32','HEAPF64','UTF8ToString']"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sSTACK_SIZE=1048576" -fexceptions)
  target_compile_options(qc_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fexceptions>)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/dist")
  foreach(asset index.html app.mjs worker.mjs molecule.mjs basis.mjs engine.mjs rhf.mjs uhf.mjs mp2.mjs electronic.mjs scan.mjs viewer.mjs plot.mjs styles.css)
    configure_file("web/${asset}" "${CMAKE_BINARY_DIR}/dist/${asset}" COPYONLY)
  endforeach()
  file(COPY web/vendor DESTINATION "${CMAKE_BINARY_DIR}/dist")
else()
  add_library(qc_core SHARED src/overlap.c src/integrals.c src/rhf.cpp src/uhf.cpp src/mp2.cpp)
endif()

target_include_directories(qc_core PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/vendor/libcint/include"
  "${CMAKE_CURRENT_BINARY_DIR}/vendor/libcint/include")
target_include_directories(qc_core SYSTEM PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen")
target_compile_definitions(qc_core PRIVATE EIGEN_MPL2_ONLY)
target_link_libraries(qc_core PRIVATE cint m)

