from setuptools import setup

setup(
    name = 'remap-gff3',
    version = '1.0',
    description = '',
    author = 'Li-Mei Chiang',
    author_email = 'dytk2134@gmail.com',
    license = 'MIT',
    install_requires = [
        'numpy',
        'gff3==0.3.0',
        'bx-python==0.7.3',
        'CrossMap'
    ],
    dependency_links = [
        'git+ssh://git@github.com:NAL-i5K/GFF3toolkit.git#egg=GFFtoolkit-1.0'
    ]
)
