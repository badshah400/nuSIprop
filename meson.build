project(
  'nuSIpropKK',
  'cpp',
  version         : '0.0.1',
  default_options : [
    'cpp_std=c++20',
    'warning_level=3',
    'strip=true',
    'b_lto=true',
    'b_pie=true',
    ]
  )

cxx = meson.get_compiler('cpp')

gomp_lib = dependency('openmp')

add_project_arguments(
  language : 'cpp'
  )

add_project_arguments('-march=native', language: 'cpp')

# if get_option('debugoutput')
#   add_project_arguments('-DDEBUG', language: 'cpp')
# endif
#
# if get_option('analyzer')
#   add_project_arguments('-fanalyzer', language: 'cpp')
# endif

incdir = include_directories('.', 'model', 'indicators')

src  = [
  'model/model.cpp',
  ]

deps = [
  dependency('fmt'),
  dependency('gsl'),
  gomp_lib,
  ]

# if get_option('profile')
#   tcmalloc_lib  = cxx.find_library('tcmalloc')
#   deps         += [ tcmalloc_lib ]
# endif

subdir('events')
subdir('chisq')

main_exe  = executable(meson.project_name(),
  src + [ 'nuflux.cpp' ],
  include_directories: incdir,
  dependencies : deps,
  install      : true,
  install_dir  : meson.project_source_root()
  )

