
# meson.build (top-level)

project('libhankel', ['c', 'cpp'],
  version: '0.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=2']
)

inc_root = include_directories('include')

# UNITY

unity_inc = include_directories('submodules/unity/src')

unity_lib = static_library(
  'unity',
  'submodules/unity/src/unity.c',
  include_directories : unity_inc,
  c_args : ['-DUNITY_INCLUDE_DOUBLE']
)

unity_dep = declare_dependency(
  link_with : unity_lib,
  include_directories : unity_inc
)

subdir('src')
subdir('include')

hankel_dep = declare_dependency(
  link_with: libhankel,
  include_directories: inc_root,
)

subdir('tests')

# rm -rf build
# meson setup build 
# meson compile -C build