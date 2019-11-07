"""
NetSurfP2

"""

from setuptools import setup


setup(
    name='netsurfp2',
    version='2.0',
    url=None,
    license='Proprietary',
    author='Michael Schantz Klausen & Martin Closter Jespersen',
    author_email='',
    description='Bioinformatic tool to predict structural protein features',
    #long_description=__doc__,
    packages=['netsurfp2', ],
    #zip_safe=False,
    platforms='any',
    install_requires=[
        'numpy>=1.15',
        'tensorflow>=1.4',
    ],
    entry_points = {'console_scripts': [
       'netsurfp2 = netsurfp2.__main__:entry',
    ]},
    include_package_data=True,
    package_data={'netsurfp2': ['model_data/*.*']},
)
