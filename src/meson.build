inc_root = include_directories('..', '../include')

# ---- Dependencies ----

# Boost.Math is used header-only, so an include path is all the build needs;
# dependency('boost') is no good here because it fails unless it also finds a
# compiled boost library, which libboost-dev does not ship. The headers are on
# the default search path on Linux, but not on macOS, where Homebrew installs
# under /opt/homebrew (Apple Silicon) or /usr/local (Intel). For anywhere else,
# pass CPPFLAGS=-I/path/to/boost.
boost_args = host_machine.system() == 'darwin' ? [
  '-I/opt/homebrew/include',
  '-I/usr/local/include',
] : []
boost_dep = declare_dependency(compile_args : boost_args)

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

