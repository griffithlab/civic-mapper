from setuptools import setup

setup(name = 'civic_mapper',
	version = '0.0.1',
	description = 'Python tool to map variants to CIViC variant interpretations',
	install_requires = [
		'requests',
		'flask',
		'transvar',
	],
	entry_points = {
		'console_scripts': ['civic-mapper = civic_mapper.command_line:main']
	},
	classifiers = [
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7',

	],
	keywords = 'civic mapper',
	url = 'http://github.com/griffithlab/civic-mapper',
	author = 'Cody Ramirez',
	author_email = 'cramirez@genome.wustl.edu',
	license = 'MIT',
	packages = ['civic_mapper'],
	include_package_data = True,
	zip_safe = False)