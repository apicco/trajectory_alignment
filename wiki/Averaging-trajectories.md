---
layout: page
title: Averaging trajectories
wikiPageName: Averaging-trajectories
menu: wiki
---

Trajectories that describe the spatio-temporal properties of an event can be aligned in space and in time together to estimate what is their average behaviour. 
The trajectories must obey to some important prerequisites, which are essential:
* The event that the trajecotories describe must be homogeneous. Trajectories describing heterogeneous events cannot be obviously averaged together.
* All the trajectories must share the same spatial and temporal units. It is not possible to align trajectories whose frames are spaced by different time intervals.
* Trajectories need to describe an event that moves in a clear direction. Unmotile events cannot be reliably aligned in space. 

The functions necessary to average the different trajectories together are defined in [average.py](https://github.com/apicco/trajectory_alignment/tree/master/trajalign) and can be imported as

	from trajalign.average import load_directory
	from trajalign.average import average_trajectories
    
**`load_directory( path , pattern = '.txt' , sep = None , comment_char = '#' , dt = None , t_unit = '' , attrs...)`** loads all the trajectories listed in `path`, which have the same `pattern`, into a list of trajectories. The trajectories in `path` must be text files with the data organised in columns separated by `sep`. When `sep` is set to _None_ the separation can be an indefinite number of white spaces. Comment lines in the trajectory files must be identified by the the `comment_char` at their beginning. `comment_char` can be a string of any length. `attrs` is used to assigne each column of the trajectory file to the right attribute and it is used to add annotations, if necessary. The [attributes](The-trajectory-class#trajectory-attributes) of a trajectory are _frames_, _t_, _coord_, _f_, _mol_, _n_, _t_err_, _coord_err_, _f_err_ and _mol_err_. 

A trajectory file can look as such:

	# This is an example of a trajectory. 
	# The file name is example.txt
	# frames x_coord y_coord fluorescence
	0 34.037666 47.123764 3.577058
	1 34.128681 47.177021 3.506692
	2 34.187073 47.090092 3.378920
	3 34.162273 47.169624 3.616615
	4 34.225525 47.163013 3.667110

The trajectory _example.txt_ can be loaded with:

	load_directory( 
		path = 'wherever the trajectory is',
		frames = 0 ,
		coord = (1 , 2) ,
		f = 3 ,
		note = 'this is an example of how to load a trajectory' ,
		experiment = 'an important one' ,
		coord_unit = 'pixels of our EM-CCD camera')

If the `dt` is given ( as a float or int ) and no _t_ attribute is later specified, then the time `.t()` attribute is created and it is computed from the frame numbers as t = ( frame_number - 1 ) * dt. As for the coordinates, a _t_unit_ attribute must be added, specifying the unit of time in which `dt` is expressed.

**`average_trajectories( trajectory_list , max_frame = 500 , output_file = 'average' , median = False )`** align all the trajectories in the list `trajectory_list` together and computes their average. 

`max_frame` is the largest frame number to be expected for the trajectories and it is used to estimate whether a trajectory is truncated at its end. Trajectories that are truncated at their start are recognised by the frame index that starts as 0. If the trajectory are complete the frame indexes should be included between 0 and `max_frame`. 
Trajectories that start with 0 are not considered when computing the start of the average trajectory, which is the average start of non-trucated trajectories.
Similarly, trajectories that end with 'max_frame' are not considered when computing the end of the average trajectory.

All trajectories are used as a reference to align all the remaining trajectories together. The best alignment, which is the one that minimise the dispersion of the aligned trajectories, is saved in a file named as in `output_file`. 

By default `average_trajectories` computes the average of the trajectories and their standard error of the mean. However, it can be useful to compute the median instead of the average, in particular when trajectories are noisy and few can be badly aligned affecting the mean. To compute the median, change the default `median = False` to `median = True`. Together with the median, the softward will compute the standard error derived from the median absolute deviation (MAD) adjusted for asymptotically normal consistency. If the standard deviation or the MAD are needed, it is sufficient to multiply the errors by the square root of the attribute `n`.

`average_trajectories` returns two trajectories, the best average and the worst average, and all the raw trajectories aligned together to compute the best average, which are oriented in space and time as the best average for comparison. 
The best average is saved as a txt file named from `output_file`. All the raw trajectories are saved in a folder names as `output_file`. The average trajectories inherit the annotations of the trajectory that is used as a reference. In addition, a new annotation _reference_file_ is added to record which is the file has been used as reference.
