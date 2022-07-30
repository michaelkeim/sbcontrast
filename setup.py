import setuptools

with open('README.txt', 'r') as file:
    long_description = file.read()

setuptools.setup(
	name='sbcontrast',
	packages=['sbcontrast'],
	version='0.1.0',
	license='MIT',
	description='Surface Brightness Limit Calculator',
	long_description=long_description,
	long_description_content_type="text/plain",
	author='Keim, M. A., van Dokkum, P., Li, J.',
	author_email='michael.keim@yale.edu',
	entry_points={"console_scripts": ["sbcontrast = sbcontrast.main:sbc"]},
	url='https://github.com/michaelkeim/sbcontrast/', 
	install_requires=['numpy', 'astropy'],
	classifiers=[
		'Development Status :: 1 - Planning',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: MIT License',
	],
	
	download_url="https://github.com/michaelkeim/sbcontrast/archive/master.tar.gz",
)
