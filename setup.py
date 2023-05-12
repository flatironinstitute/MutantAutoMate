from setuptools import setup, find_packages

setup(
    name='packaging_demo',
    version='1.0',
    author='Akshada Pradhan',
    author_email='asameerpradhan@flatironinstitute.org',
    description='A demo package',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
    ],
)
