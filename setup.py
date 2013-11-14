from distutils.core import setup

setup(
    name='GroopM',
    version='0.2.10.12',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['groopm', 'groopm.test'],
    scripts=['bin/groopm'],
	url='http://pypi.python.org/pypi/GroopM/',
    license='LICENSE.txt',
    description='Metagenomic binning suite',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.6.1",
        "scipy >= 0.10.1",
        "matplotlib >= 1.1.0",
        "tables >= 2.3",
        "pysam >= 0.6",
        "PIL >= 1.1.7",
        "BamTyper >= 0.2.6"
    ],
)
