from setuptools import setup

setup(name='oatool',
        version='0.0.1',
        description='Online analysis toolkit in Python',
        author='Benjamin Naecker',
        author_email='bnaecker@fastmail.com',
        url='https://github.com/baccuslab/oatool',
        long_description='''
            This package provides a command-line tool for running basic
            online analyses of retinal data during an experiment. It currently
            allows experimenters to easily compute online receptive fields,
            both via reverse-correlation for intracellular recordings and via
            spike-triggered averaging for extracellular. It also provides an
            extensible library (oalib) for writing new analyses as needed.
            ''',
        classifiers=[
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering',
            'License :: OSI Approved :: GPL',
            'Operating System :: OS Indepedent',
            'Programming Language :: Python :: 3'
        ],
        install_requires=[
            'pyret>=0.5.4',
            'bldsclient>=0.0.1'
        ]
    )
