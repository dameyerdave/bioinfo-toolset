import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='bioinfo_toolset',
    version=1.4,
    author='David Meyer',
    author_email='dameyerdave@gmail.com',
    description='Python Rules Evaluator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.5',
    url='https://github.com/dameyerdave/bioinfo-toolset',
    packages=setuptools.find_packages(exclude='test'),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'liftover',
        'click',
        'friendlylog',
        'termcolor',
        'tqdm',
        'hgvs',
        'jsonschema',
        'docker',
        'ga4gh.vrs[extras]'
    ],
    entry_points={
        'console_scripts': {
            'bit = bioinfo_toolset.bit:cli'
        }
    }
)
