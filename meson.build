
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
  'external_libs_based_work/qwe/tanhsinh.c',
  'src/utils/analytical_form_factors.c'
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



# Export for tests/consumers:
libhankel_dep = declare_dependency(
  link_with: libhankel,
  include_directories: inc,
  dependencies: [gsl_dep],     # <-- carries -lgsl -lgslcblas -lm
)


tests = [
  'tests/test_DHT_hankel.c', 
  'tests/test_analytical.c',
]

foreach t : tests
  exe = executable(
    t.replace('.c', '').replace('tests/', ''),
    t,
    dependencies : [libhankel_dep, unity_dep]
  )
  test(t, exe)
endforeach



# run with meson test -C build
#test('mylib works', test_exe, args: [meson.current_source_dir() / 'tests/data.csv'])