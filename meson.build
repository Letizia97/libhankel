# meson.build (top-level)

project(
  'libhankel', 
  ['c', 'cpp'],
  meson_version : '>=1.4.0',
  version: '0.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=2']
)

inc_root = include_directories('include')

cc = meson.get_compiler('c')

# Apple clang rejects a call to an undeclared function outright while GCC only
# warns, so make it an error on both: otherwise a missing declaration builds
# cleanly on Linux and only breaks once it reaches macOS.
add_project_arguments(
  cc.get_supported_arguments('-Werror=implicit-function-declaration'),
  language : 'c',
)

subdir('src')
subdir('include')

hankel_dep = declare_dependency(
  link_with: libhankel,
  include_directories: inc_root,
)

if get_option('tests')
  subdir('tests')
endif



# TO BUILD:

# rm -rf build

# Build library only
# meson setup build

# OR Build library + tests
# meson setup build -Dtests=true

# meson compile -C build