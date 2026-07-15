inc_root = include_directories('..', '../include')

# ---- Main C library ----
sources = files(
  'hankel_DHT.c', 
  'form_factors.c', 
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
  dependencies: [],          
  install: true,
)

