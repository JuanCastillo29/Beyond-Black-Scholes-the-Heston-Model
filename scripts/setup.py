# -*- coding: utf-8 -*-
"""
Created on Wed Aug 01 10:35:41 2025

@author: juanc
"""
from distutils.core import setup, Extension
import os
import shutil

# Configura las extensiones con optimización
lib_Heston = Extension(
    name='Heston_lib',
    sources=['Heston_MonteCarlo_simulation.c'],
    extra_compile_args=['/O2'],  # Optimización para velocidad
)



lib_HestonFFT = Extension(
    name='HestonFFT_lib',
    sources=['Heston_MonteCarlo_FFT.c'],
    extra_compile_args=['/O2'],  # Optimización para velocidad
)

# Crear un setup temporal para compilar ambas librerías
setup(
    name='custom_dlls',
    version='1.0',
    description='Compila librerías Heston',
    ext_modules=[lib_Heston, lib_HestonFFT]
)

# Mueve los .dll a la raíz del proyecto
def mover_dlls():
    build_dir = os.path.join(os.getcwd(), 'build')
    for root, dirs, files in os.walk(build_dir):
        for file in files:
            if file.endswith('.dll'):
                src = os.path.join(root, file)
                dst = os.path.join(os.getcwd(), file)
                shutil.copyfile(src, dst)
                print(f"{file} compilado y movido a {dst}")

if __name__ == "__main__":
    mover_dlls()