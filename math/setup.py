import pybind11
from distutils.core import setup, Extension

ext_modules = [
    Extension(
        'numera',
        ['math.cpp', 'parser.cpp'],
        include_dirs=['D:\pybind11-master\install\include'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
]

setup(
    name='numera',
    version='0.0.1',
    author='TheCakeMan1',
    author_email='lilo.lilo565@gmail.com',
    #description='pybind11 extension',
    ext_modules=ext_modules,
    requires=['pybind11']
)
