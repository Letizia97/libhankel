inc_root = include_directories('..', '../include')

cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')

# ---- Dependencies ----

# The maths library is a separate object on Linux but part of libSystem on
# macOS, so it is looked up rather than hard-coded as '-lm'.
m_dep = cc.find_library('m', required : false)

# Only the header-only Boost.Math special functions are used.  Debian and
# Ubuntu put those headers on the default include path, so a plain compile
# check finds them; Homebrew keeps them under its own prefix, which the
# compiler does not search, so ask brew where they are before giving up.
boost_header = 'boost/math/special_functions/bessel.hpp'
boost_dep = dependency('boost', required : false)
boost_inc = []
boost_found = boost_dep.found() or cxx.has_header(boost_header)

if not boost_found
  brew = find_program('brew', required : false)
  if brew.found()
    brew_boost = run_command(brew, '--prefix', 'boost', check : false)
    if brew_boost.returncode() == 0
      candidate = brew_boost.stdout().strip() / 'include'
      if cxx.has_header(boost_header, args : ['-I' + candidate])
        boost_inc = include_directories(candidate)
        boost_found = true
      endif
    endif
  endif
endif

if not boost_found
  error(
    'Boost.Math headers not found. Install Boost - "apt install libboost-dev" on ' +
    'Debian/Ubuntu, "brew install boost" on macOS - or point the compiler at it ' +
    'with BOOST_ROOT or CPPFLAGS.'
  )
endif

lib_deps = []
if boost_dep.found()
  lib_deps += boost_dep
endif
if m_dep.found()
  lib_deps += m_dep
endif

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
  include_directories: [inc_root, boost_inc],
  dependencies: lib_deps,
  install: true,
)

