
test_utils_lib = static_library(
  'test_utils',
  ['test_utils.c'],
  include_directories: include_directories('.'),
  dependencies: [unity_dep] 
)

test_utils_dep = declare_dependency(
  link_with: test_utils_lib,
  include_directories: include_directories('.')
)
