gsl_dep = dependency('gsl', required: true, method: 'pkg-config')
inc_root = include_directories('..', '../include')
ogata_inc = include_directories('../external_libs_based_work/Ogata')


# ---- Ogata C++ sources ----
ogata_cpp_src = files(
  # List the actual Ogata .cpp files you want to compile
  '../external_libs_based_work/Ogata/FTB.cpp',
)

# Compile them as a static or object library (object library avoids archive)
ogata_objlib = static_library(
  'ogata_cpp',
  sources: ogata_cpp_src,
  include_directories: [ogata_inc],
)

ogata_dep = declare_dependency(
  link_with: ogata_objlib,
  include_directories: [ogata_inc],
)

# ---- Your C++ wrapper that exposes a C ABI ----
wrapper_cpp = files('wrappers/wrapper_FBT.cpp')

wrapper_lib = static_library(
  'FBT',
  sources: wrapper_cpp,
  include_directories: [inc_root, ogata_inc],
  dependencies: [ogata_dep, gsl_dep],
  install: false,
)

wrapper_dep = declare_dependency(
  link_with: wrapper_lib,
  include_directories: [inc_root],
)


# ---- Main C library ----
sources = files(
  'hankel_DHT.c', 
  'form_factors.c', 
  'hankel_FBT_Ogata.c', 
  'hankel_DE_quadrature.c', 
  'hankel_QWE.c', 
  '../external_libs_based_work/DE-quadrature/intde.c', 
  '../external_libs_based_work/qwe/qwe_Chave.c', 
  '../external_libs_based_work/qwe/qwe_Key.c', 
  'utils/sasfit_integrate.c',
  '../external_libs_based_work/qwe/tanhsinh.c',
  'utils/analytical_form_factors.c',
)

libhankel = library(
  'libhankel',                             # builds liblibhankel.so/.a
  sources: sources,
  include_directories: [inc_root],
  dependencies: [wrapper_dep, gsl_dep],          
  install: true,
)

