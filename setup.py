from setuptools import setup, find_packages
#from setuptools.command.develop import develop
#from setuptools.command.install import install

setup(
        name="PyLumerical",
        version="0.1.0",
        packages=find_packages(exclude=['*test']),
        author="Gareth Sion Jones",
        author_email="gareth.jones@materials.ox.ac.uk",
        python_requires='>3.4.0',
        install_requires=['numpy', 'matplotlib']
    ) 
