from setuptools import setup, Extension

module_numera = Extension('numera', sources=['main.c', 'harmonic_c.c', 'struct.c'])
#module_harmonic = Extension('numera.harmonic', sources=['harmonic_c_main.c', 'harmonic_c.c', 'struct.c'])

setup(
    name='numera',
    version='1.0',
    description='Discrete Fourier Transform Module',
    ext_modules=[module_numera]
)


'''import pybind11
from distutils.core import setup, Extension

ext_modules = [
    Extension(
        'numera',
        ['math.cpp', 'parser.cpp', 'harmonic.cpp', 'vector.cpp', 'val.cpp', 'vector_c.c'],
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
)'''