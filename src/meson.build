inc_root = include_directories('..', '../include')

# ---- Dependencies ----

# Only the header-only Boost.Math special functions are used, so an include path
# is all that is needed. dependency('boost') is no good here: it fails unless it
# also finds a compiled boost library, which headers-only installs (libboost-dev)
# do not ship. The headers are on the default search path on Linux; on macOS
# Homebrew puts them under /opt/homebrew (Apple Silicon) or /usr/local (Intel).
cxx = meson.get_compiler('cpp')
boost_header = 'boost/math/special_functions/bessel.hpp'

boost_args = []
boost_found = cxx.has_header(boost_header)
foreach prefix : ['/opt/homebrew/include', '/usr/local/include']
  if not boost_found and cxx.has_header(boost_header, args : ['-I' + prefix])
    boost_args = ['-I' + prefix]
    boost_found = true
  endif
endforeach

if not boost_found
  error(
    'Boost.Math headers not found. Install them with "sudo apt install ' +
    'libboost-dev" or "brew install boost", or point the build at an existing ' +
    'copy with CPPFLAGS=-I/path/to/boost.',
  )
endif

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

