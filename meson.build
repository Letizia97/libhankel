# meson.build (top-level)

project(
  'libhankel', 
  ['c', 'cpp'],
  meson_version : '>=1.4.0',
  version: '0.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=2']
)

inc_root = include_directories('include')

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