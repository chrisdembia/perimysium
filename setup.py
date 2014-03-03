from setuptools import setup

setup(
        name='perimysium',
        author='Chris Dembia',
        author_email='chris530d@gmail.com',
        version='0.1',
        description='Client code code for OpenSim.',
        extras_require={'doc': ['sphinx', 'numpydoc']},
        packages=['perimysium'],
        scripts=['bin/stoplot'],
        url='http://github.com/fitze/perimysium',
        license='LICENSE.txt',
        long_description=open('README.txt').read(),
        )

