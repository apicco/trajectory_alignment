---
layout: page
title: Installation
wikiPageName: Installation
menu: wiki
---

`trajalign` uses exclusively [Python 3.X](https://docs.python.org/3/). Before proceeding with the installation make sure to have following libraries:

* [numpy](http://www.numpy.org/) fundamental library for scientific computing;
* [sklearn](http://scikit-learn.org/stable/) Scikit machine learning library.

On a mac: if those libraries are not installed, you can install them with [macport](https://www.macports.org/):
	
	sudo port install py35-numpy py35-scikit-learn #or the latest versions you will find..
 
`trajalign` uses also the following libraries that should be present as defaults in normal python installations:

* [math](https://docs.python.org/2/library/math.html),
* [os](https://docs.python.org/2/library/os.html),
* [copy](https://docs.python.org/2/library/copy.html),
* [warnings](https://docs.python.org/2/library/warnings.html).

To plot the trajectories of the examples we use [matplotlib](http://matplotlib.org/). This library is not necessary to use the functionalities of the package.

To install the `trajalign` package, please download the latest version of the package from [here](https://github.com/apicco/trajectory_alignment/archive/master.zip) or visit the [github repository](https://github.com/apicco/trajectory_alignment).

In the folder that you downloaded you will find a `setup.py` file. From the shell type:

	sudo python3.5 setup.py install

The `trajalign` package is now ready to use. Refer to the [documentation](Home) and to these [examples](Examples), which are also present in the folder you downloaded, to learn how to use the package.

Have fun!


