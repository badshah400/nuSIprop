# project(
#   'KKEvents',
#   'cpp',
#   default_options : [
#     'cpp_std=c++20',
#     'warning_level=3',
#     'strip=true',
#     'b_lto=true',
#     'b_pie=true',
#     ]
# )

# cxx = meson.get_compiler('cpp')
event_deps = [
  dependency('tbb'),
  dependency('threads'),
  dependency('fmt')
]

event_exe  = executable('chisq.x',
  sources: [ 'chisq.cpp' ],
  include_directories: include_directories('.', '../events/'),
  dependencies : event_deps,
  install      : true,
  install_dir  : meson.project_source_root()
  )

