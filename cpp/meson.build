project('rd2d', 'cpp')

#  default_options : ['cpp_std=c++11'])

boost_dep = dependency('boost', modules : ['timer', 'system', 'program_options'])

executable('rd2d', 'rd2d.cpp', dependencies : boost_dep, cpp_args : ['--std=c++11', '-march=native', '-O3'])
