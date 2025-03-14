import setuptools

setuptools.setup(
    name='SpectraAnalyser',
    version='1.1',    
    description='Code to analyse spectral data.',
    url='https://github.com/AllisonPessoa/luminescence-data-analyser',
    author='Allison Pessoa',
    author_email='allison.pessoa@ufpe.br',
    license='GNU GPL',
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    install_requires=['matplotlib',
                      'scipy', 'numpy', 'PrettyTable', 'uncertainties'
                      ],
    package_data={'pcf_lib':['*',]}
)
