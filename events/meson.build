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

event_exe  = executable('events.x',
  sources: [ 'flux_integration.cpp', 'events.cpp' ],
  include_directories: include_directories('.'),
  dependencies : event_deps,
  install      : true,
  install_dir  : meson.project_source_root()
  )

