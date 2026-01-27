
tests = [
  'test_analytical.c', 
  'test_hankel_DHT.c', 
  'test_hankel_QWE.c'
]
test_inc = include_directories('.')   


utils_lib = static_library(
  'utils',
  'utils.c',
  dependencies : [unity_dep],
  include_directories: test_inc,
)

foreach t : tests
  exe = executable(
    t.replace('.c', '').replace('tests/', ''),
    t,
    dependencies : [hankel_dep, unity_dep],
    include_directories: test_inc,
    link_with: utils_lib,
  )
  test(t, exe)
endforeach


# run with meson test -C build
#test('mylib works', test_exe, args: [meson.current_source_dir() / 'tests/data.csv'])
