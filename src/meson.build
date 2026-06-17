gsl_dep = dependency('gsl', required: true, method: 'pkg-config')
inc_root = include_directories('..', '../include')
ogata_inc = include_directories('../external_libs/Ogata')


# ---- Ogata C++ sources ----
ogata_cpp_src = files(
  # List the actual Ogata .cpp files you want to compile
  '../external_libs/Ogata/FTB.cpp',
)

# Compile as a static or object library 
ogata_objlib = static_library(
  'ogata_cpp',
  sources: ogata_cpp_src,
  include_directories: [ogata_inc],
)

ogata_dep = declare_dependency(
  link_with: ogata_objlib,
  include_directories: [ogata_inc],
)

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
  'hankel_transform.c',
  '../external_libs/DE-quadrature/intde.c', 
  '../external_libs/qwe_Chave/qwe_Chave.c', 
  '../external_libs/qwe_Key/qwe_Key.c', 
  'utils/sasfit_integrate.c',
  '../external_libs/utils/tanhsinh.c',
  'utils/analytical_form_factors.c',
  'utils/sf_functions.c',
  'utils/boost_bessel.cpp'
)

libhankel = shared_library(
  'libhankel',                           
  sources: sources,
  include_directories: [inc_root],
  dependencies: [wrapper_dep, gsl_dep],          
  install: true,
)

