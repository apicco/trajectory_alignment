---
layout: page
title: The trajectory class
wikiPageName: The-trajectory-class
menu: wiki
---

Here, all the elementss of the trajectory class are listed. The trajectory class is defined in [format/traj.py](https://github.com/apicco/trajectory_alignment/tree/master/format) and can be imported as

	from trajalign.traj import Traj

## Trajectory definition

**`Traj(annotations)`** creates a new empty trajectory. Annotations are an optional dictionary allowing whatever entries are needed about the trajectory (name, date, experiment,..)

`t = Traj( what = 'my first trajectory' , folder = 'home/foo/' , experiment = 'very important' , mood = 'today is a beautiful day' )`

Annotations can be read and modified at anytime (see `.annotations()`).

`print(t)` prints the trajectory content and its annotations.

## Trajectory attributes

A trajectory contains several attributes:

* **`.frames()`** is an array of frame indexes of the original image. If the time is present, then `.frames()` and `.t()` must have the same length.

* **`.t()`** is an arrey of time values associated to each frame, if any. The array of `.t()` must have the same length as `.frames()`.

* **`.coord()`** is a matrix of the spatial coordinate of the fluorescent spot centroid.

* **`.f()`** an array of fluorescence intensities measured from the spots.

* **`.mol()`** an array of number of molecules, if known.

* **`.n()`** if traj() is an average trajectory, then `.n()` outputs an array containing the number of trajectories used to compute the average at each datapoint.

`.t()`, `.coord()`, `.f()` and `.n()` have related attributes containing the measured errors, if known. The errors are called with `.t_err()`, `.coord_err()`, `.f_err()` and `.mol_err()` respectively.

**`.attributes()`** lists the attributes of the trajectory that are not empty.

**`.input_values( attribute , array , unit = '' )`** sets attribute values.

`d.input_values( 't'  ,  [ 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1 ]  ,  unit = 's' )`

`d.input_values( 'coord'  ,  [[ 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1] , [1 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ]] , unit = 'pixels' )`

**`.norm_f()`** normalise the fluorescence intensity `.f()` between 0 and 1.

## Properties of the trajectory class:

Each attribute of the trajectory allows integer arguments to index elements in its array, which are then returned.

`print( d.t( 0 , 1 , 2 ) )`

`print( d.t( [ 0 , 1 , 2 ] ) )`

`print( d.t( range( 0 , 3 ) ) )`

or (deprecated) `print( d.t()[ range( 0 , 3 ) ] )`

`print( d.coord( 0 , 1 , 2 ) )`

`print( d.coord( [ 0 , 1 , 2 ] ) )`

`print( d.coord( range( 0 , 3 ) ) )`

or (deprecated) `print( [ d.coord()[ 0 ][ range( 0 , 3 ) ] , d.coord()[ i ][ range( 0 , 3 ) ] ] )`

**`.annotations( annotation = None , string = '')`**: output the dictionary of the annotations associate to the trajectory (equivalent to `.__dict__()`). If an annotation is inputed it changes the value of the annotation with string. If the annotation is not defined it defines it.

Annotations can be outputted:

`print( d.annotations()[ 'what' ])`

They can be added

`d.annotations( 'goals' , 'Understand the meaning of live' )`

or modified

`d.annotations( 'mood' , 'all this seems pretty complicated. But it is not!' )`

which is equivalent to (deprecated) `d.annotations()[ 'mood' ]= 'all this seems pretty complicated. But it is not!'`

**`.start( t = None )`** returns the start of the trajectory if the time attribute is defined. If t is specified then it returns the trajectory attributes starting from t.

`print( d.start() )`

`print( d.start( t = 0.5 ) )`

**`.end( t = None )`** returns the end of the trajectory if the time attribute is defined. If t is specified then it returns the trajectory attributes ending at t.

`print( d.end() )`

`print( d.end( t = 0.5 ) )`

**`.len( d )`** outputs the length of the trajectory (i.e. the length of the arrays in the trajectory attributes).

**`.extract( items )`** extracts the listed items from the trajectory attributes.

`print( d.extract( 2 , 3 , 8 )`

`print( d.extract( [ 2 , 3 , 8 ] )`

`print( d.extract( range( 0 , len( d ) ) ) )` is the identity to `d`.


## Operations on the trajectories

**`.fill()`** fills attributes with NaN where there are missing frames.

**`.rotate( angle , angle_err = 0 )`** rotates the coordinates by an angle expressed in radiants and propagate the errors accordingly if `.coord_err()` is not empty or if angle_err is not 0.

**`.translate( v , v_err = ( 0 , 0 ) )`** translate the trajectory by a vector v. If v_err is not 0 it propagates the error accordingly.

**`.lag( time_shift )`** shift in time of the trajectory by time_shift, in the trajectory time units.

**`.center_of_mass()`** translate the trajectory so that its cententer of mass sits on the origin.

**`.lifetime( round = 2 )`** computes the lifetime of the trajectory and rounds it to 'round' number of float digits.

**`.time( delta_t , unit )`** times the trajectory using the frame information. Frames are spaced by delta_t of given unit. `.time()` operates only if the `.t()` attribute is empty. 

## Input and output

**`.save( filename )`** saves the trajectory as txt to the filename.

`t.save( filename = 'my_first_trajectory.txt' )`

**`.load( filename , sep = None , comment_char = '#' , attribute_names... )`**: loads data from a txt table. Data must be ordered in columns. Columns can be separated by spaces or tabs or comas (for .csv files). The separator can be entered in sep as a string. The default for sep is None, which will recognise columns separated by an arbitrary number of spaces. The strings starting with the comment_char string will be disregarded. Attribute_names associate each column to the right attribute. Attribute_names can be only: 'frames', 't', 'coord', 'f', 'n', 'mol','t_err', 'coord_err', 'f_err', 'mol_err'.

	t = Traj() #empty trajectory
	t.load( 'my_first_trajectory.txt' , frames = 0 , t = 1 , coord = ( 2 , 3 ) )

Note that 'coord' requires two values and the column indexing starts from 0.
