from setuptools import setup

setup(
        name='epimysium',
        author='Chris Dembia',
        author_email='chris530d@gmail.com',
        version='0.1',
        description='Client code code for OpenSim.',
        extras_require={'doc': ['sphinx', 'numpydoc']},
        packages=['epimysium'],
        url='http://github.com/fitze/epimysium',
        license='LICENSE.txt',
        long_description=open('README.txt').read(),
        )

