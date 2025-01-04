from setuptools import setup, find_packages

setup(
    name="rumaps",
    version="0.1.1",
    description='Refined UMAPS (RUMAPS) is a Python package crafted for generating polished UMAP visualizations, offering streamlined labeling and interactive exploration specifically for scanpy adata object.',
    author="Sai Rama Sridatta Prakki",
    author_email="ramadatta.88@gmail.com",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/ramadatta/rumaps",
    packages=find_packages(),
    install_requires=[
        "scanpy",
        "datamapplot",
        "matplotlib",
        "numpy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
