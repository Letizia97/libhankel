
subdir('utils')

tests = [
  'test_analytical.c', 
  'test_hankel_DHT.c', 
  'test_hankel_QWE_Chave.c',
  'test_hankel_QWE_Key.c',
  'test_integrate_ctm.c',
  'test_hankel_transform.c',
  'test_form_factors.c',
  'test_hankel_DE_Ogata.c',
  'test_hankel_DE_Ooura.c'
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
