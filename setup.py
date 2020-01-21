import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dna_designer2", # Replace with your own username
    version="0.0.2",
    author="Keoni Gandall",
    author_email="koeng101@gmail.com",
    description="A DNA design package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://github.com/koeng101/dna_designer2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

