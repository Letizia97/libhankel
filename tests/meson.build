
tests = ['test_analytical.c', 'test_DHT_hankel.c']
test_inc = include_directories('.')   


test_common_lib = static_library(
  'test_common',
  'test_common.c',
  dependencies : [unity_dep],
  include_directories: test_inc,
)

foreach t : tests
  exe = executable(
    t.replace('.c', '').replace('tests/', ''),
    t,
    dependencies : [hankel_dep, unity_dep],
    include_directories: test_inc,
    link_with: test_common_lib,
  )
  test(t, exe)
endforeach


# run with meson test -C build
#test('mylib works', test_exe, args: [meson.current_source_dir() / 'tests/data.csv'])