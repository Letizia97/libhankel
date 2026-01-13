
# meson.build (top-level)

project('libhankel', ['c', 'cpp'],
  version: '0.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=2']
)

gsl_dep = dependency('gsl', required: true, method: 'pkg-config')

# Include and library directories for your local external_libs_based_work
inc = include_directories('include')
ogata_inc = include_directories('external_libs_based_work/Ogata')


# ---- Ogata C++ sources ----
ogata_cpp_src = files(
  # List the actual Ogata .cpp files you want to compile
  'external_libs_based_work/Ogata/FTB.cpp',
)

# Compile them as a static or object library (object library avoids archive)
ogata_objlib = static_library(
  'ogata_cpp',
  sources: ogata_cpp_src,
  include_directories: [ogata_inc],
  # You can add cpp_args if the vendor requires specific flags
)

ogata_dep = declare_dependency(
  link_with: ogata_objlib,
  include_directories: [ogata_inc],
)

# ---- Your C++ wrapper that exposes a C ABI ----
wrapper_cpp = files('src/wrappers/wrapper_FBT.cpp')

wrapper_lib = static_library(
  'FBT',
  sources: wrapper_cpp,
  include_directories: [inc, ogata_inc],
  dependencies: [ogata_dep, gsl_dep],
  install: false,
)

wrapper_dep = declare_dependency(
  link_with: wrapper_lib,
  include_directories: [inc],
)

# ---- Your main C library that calls the C ABI ----
libhankel_c = files(
  'src/hankel_DHT.c', 
  'src/form_factors.c', 
  'src/hankel_FBT_Ogata.c', 
  'src/hankel_DE_quadrature.c', 
  'src/hankel_QWE.c', 
  'external_libs_based_work/DE-quadrature/intde.c', 
  'external_libs_based_work/qwe/qwe_Chave.c', 
  'external_libs_based_work/qwe/qwe_Key.c', 
  'src/utils/sasfit_integrate.c',
  'external_libs_based_work/qwe/tanhsinh.c'
)

libhankel = library(
  'libhankel',                             # builds liblibhankel.so/.a
  sources: libhankel_c,
  include_directories: [inc],
  dependencies: [wrapper_dep, gsl_dep],           # pulls in wrapper and transitively ogata_cpp
  link_language: 'cpp',                  # important: ensure C++ linker is used
  install: true,
)

install_headers('external_libs_based_work/Ogata/FBT.h', subdir: 'libhankel')



# Export for tests/consumers:
libhankel_dep = declare_dependency(
  link_with: libhankel,
  include_directories: inc,
  dependencies: [gsl_dep],     # <-- carries -lgsl -lgslcblas -lm
)

test_exe = executable(
  'test_libhankel',
  ['tests/test_libhankel.c'],
  dependencies: [libhankel_dep],
)

# run with meson test -C build
test('mylib works', test_exe)