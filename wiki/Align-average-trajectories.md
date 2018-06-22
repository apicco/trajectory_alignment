---
layout: page
title: Align average trajectories
wikiPageName: Align-average-trajectories
menu: wiki
---

Average trajectories can be aligned in space and in time together by chosing one protein, and its average trajectory, as a spatial and temporal landmark. 
For simplicity, we call this protein _reference protein_ and its trajectory _reference trajectory_, while the protein that needs to be alinged is referred to as _target proteins_ and its dynamics are described by the _target trajectory_.

The functions necessary to average the different trajectories together are defined in 
To align the target trajectory to the reference trajectory we need to import the function `align` in [align.py](https://github.com/apicco/trajectory_alignment/tree/master/trajalign) and the function `load_directories` in [average.py](https://github.com/apicco/trajectory_alignment/tree/master/trajalign):

	from trajalign.align import align
	from trajalign.average import load_directory

`load_directory` is used to load all the trajectories 
that are acquired simultaneously for both the target and the reference proteins in two distinct lists. 
These trajectories will be used to compute the transformation
that alinged the target and reference trajectories together.
It is important that the trajectories in these lists are ordered
in the same way, so that the first trajectory of the reference protein
in one list and the first trajectory of the target protein in the other 
list have been acquired in the same event. Also, it is important the  
relative position between the trajectories of the reference and target proteins
have been corrected for chromatic aberrations. See the [example](Example-of-trajectory-alignment.md) 
for details on how to load the trajectories that are used to compute
the alignment of the reference and target average trajectories.

**`align( path_target , path_reference , ch1 , ch2 , fimax1 = False , fimax2 = False )`** align the target average trajectory, identified by the path 
_path_target_ to the reference average trajectory, identified by the path _path_reference_. _ch1_ and _ch2_ are
the lists of the trajectories for the reference and target protein that have been acquired simultaneosly and which 
are used to compute the alignment. _fimax1_ and _fimax2_ are logic values that instruct align whether to use only the 
trajectory information up to the peak in fluorescence intensity (_True_) or to use the complete trajectory information (_False_, which is the default).
It is important when trajectories are particularly noisy and dynamic toward the end of their lifetime, as it is the case for Abp1 in endocytosis. 

Aligned trajectories are saved with the same name as the target trajectory (path_target), followed by "aligned". 
These trajectories store also new  annotations: the name of the target trajectory, the name of the reference trajectory and the transformation and its error that aligns the target trajectory to the reference trajectory. The transformation is computed as the median and the errors are the standard errors computed from the MAD adjusted for asymptotically normal consistency (_k = 1.4826_). 
