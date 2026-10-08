inc_root = include_directories('..', '../include')

# ---- Dependencies ----

# Boost.Math is used header-only, so an include path is all the build needs;
# dependency('boost') is no good here because it fails unless it also finds a
# compiled boost library, which libboost-dev does not ship. The headers are on
# the default search path on Linux, but not on macOS (Homebrew) or Windows
# (vcpkg). For anywhere else, pass CPPFLAGS=-I/path/to/boost.
fs = import('fs')

inc_dirs = [inc_root]
if host_machine.system() == 'darwin'
  foreach d : ['/opt/homebrew/include', '/usr/local/include']
    if fs.is_dir(d)
      inc_dirs += [include_directories(d)]
    endif
  endforeach
elif host_machine.system() == 'windows'
  inc_dirs += [include_directories('C:/vcpkg/installed/x64-windows/include')]
endif

# Separate library on Linux; part of libSystem on macOS and the C runtime on
# Windows, so required: false avoids a build failure on those platforms.
m_dep = cc.find_library('m', required: false)

# ---- Main C library ----
sources = files(
  'hankel_DHT.c',
  'form_factors.c',
  'hankel_DE_quadrature.c',
  'hankel_G0.c',
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
  'utils/boost_bessel.cpp',
  'utils/interp_cubic.cpp',
  'utils/interp_linear.c',
  'utils/interp_loglinear.c',
  'utils/tabulated_ff.c'
)

libhankel = shared_library(
  'libhankel',
  sources: sources,
  include_directories: inc_dirs,
  c_args:   ['-DLIBHANKEL_BUILDING_DLL'],
  cpp_args: ['-DLIBHANKEL_BUILDING_DLL'],
  dependencies: [m_dep],
  install: true,
)

