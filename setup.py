from setuptools import setup, find_packages

# ~ import sys
# ~ sys.path.remove('/Users/joseivan/pyflowmaps/pyflowmaps')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='pyflowmaps',
    version="0.0.2",
    author=["Jose Ivan Campos Rozo",'Santiago Vargas Dominguez'],
    author_email="hypnus1803@gmail.com",
    description='Package to infer horizontal proper motions',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Hypnus1803/pyflowmaps",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=find_packages(),#where="/Users/joseivan/opt/anaconda3/lib/python3.7/site-packages/pyflowmaps"),
    python_requires=">=3.7",
    install_requires=[
    'numpy >= 1.19','scipy >= 1.5','astropy >= 4.2','sunpy >= 3.0'
    ],
    )
