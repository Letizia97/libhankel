inc_root = include_directories('..', '../include')

# ---- Dependencies ----

# Only the header-only Boost.Math special functions are used. Meson searches
# the Homebrew prefixes as well as the system ones, so this covers Linux and
# macOS; set BOOST_ROOT for an install in an unusual place.
boost_dep = dependency('boost')

# Separate library on Linux, part of libSystem on macOS, hence the lookup
# rather than a hard-coded '-lm'.
m_dep = cc.find_library('m')

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
  'utils/validate_x.c',
  'utils/boost_bessel.cpp'
)

libhankel = shared_library(
  'libhankel',
  sources: sources,
  include_directories: [inc_root],
  dependencies: [boost_dep, m_dep],
  install: true,
)

