from setuptools import setup, find_packages

setup(
    name='phallett',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'phallett=phallett.__main__:main',  # adjust if your CLI entry point is different
        ],
    },
    install_requires=[
        # list required Python packages here if any
    ],
    author='Nathalia Portilla',
    description='A toolkit for genomic distance calculation and graphing.',
    license='MIT',
)
