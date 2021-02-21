from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
	long_description = f.read()
with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as r:
	requirements = r.read()

setup(
	name = 'pynida',
	packages = ['pynida','pynida.classes'],
	version = '0.1.16',
	license='GPLv3',
	description = 'In-situ nanoindentation data analysis',
	long_description=long_description,
	long_description_content_type='text/markdown',
	author = 'Vasily Lebedev, Anton Poluboiarinov, Daniil Kozlov, Robert Sakaev',
	author_email = 'vasily.lebedev@ul.ie',
	url = 'https://github.com/LebedevV/pynida',
	download_url = 'https://github.com/LebedevV/pynida',
	keywords = ['Nanomechanics', 'Microscopy', 'in-situ'],
	install_requires=requirements ,
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Operating System :: OS Independent',
		'Programming Language :: Python',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Topic :: Scientific/Engineering',
		'Programming Language :: Python :: 3',
		],
)
