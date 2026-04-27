
subdir('utils')

tests = [
  'test_hankel_FBT.c', 
]
test_inc = include_directories('.')   

foreach t : tests
  exe = executable(
    t.replace('.c', '').replace('tests/', ''),
    t,
    dependencies : [
      hankel_dep, 
      unity_dep, 
      test_utils_dep
    ],
    include_directories: test_inc,
  )
  test(t, exe)
endforeach


# run with meson test -C build
#test('mylib works', test_exe, args: [meson.current_source_dir() / 'tests/data.csv'])
