import setuptools
from setuptools import find_packages
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scHiCStackL",
    version="0.0.3",
    author="YingFu Wu",
    author_email="18821658087@163.com",
    description="scHiCStackL Package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    include_package_data=True,
    package_dir={"": "src"},
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        "": ["*.*"],
        # And include any *.msg files found in the "hello" package, too:
        "Files": ["Base/*.txt",
                  "PCA_file/human/626/*.txt",
                  "PCA_file/human/2655/*.txt",
                  "PCA_file/mouse/178/*.txt",
                  "test_samples/*.*",
                  "Base/*.txt"],
        "parameters":["178/*.json",
                      "626/*.json",
                      "2655/*.json"]
    },
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)