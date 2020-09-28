set(EXAMPLES
  example_amatrix_bem3d
  example_hmatrix_bem3d
  amatrix_bem3d_laplace
  amatrix_bem3d_helmholtz
  )

include_directories(${PROJECT_SOURCE_DIR}/Library)

foreach(t ${EXAMPLES})
  set(example_executable ${t})
  add_executable(${example_executable} ${example_executable}.c)
  target_link_libraries(${example_executable} H2LIB ${gcov} -lm)
endforeach()
