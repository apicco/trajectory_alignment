# All the software here is distributed under the terms of the GNU General Public License Version 3, June 2007. 
# Trajalign is a free software and comes with ABSOLUTELY NO WARRANTY.
# 
# You are welcome to redistribute the software. However, we appreciate is use of such software would result in citations of 
# Picco, A., Kaksonen, M., _Precise tracking of the dynamics of multiple proteins in endocytic events_,  Methods in Cell Biology, Vol. 139, pages 51-68 (2017)
# http://www.sciencedirect.com/science/article/pii/S0091679X16301546
# 
# Author: Andrea Picco (https://github.com/apicco)
# Year: 2017

import os 
from trajalign.traj import Traj
import copy as cp
import numpy as np
import warnings as wr

from sklearn import linear_model

def header( version = 1.90 , year = 2020 , printit = True ) :

	if printit :

		print('|-----------------------------------------------------|')
		print('| Trajalign version ' + str( version ) +' Copyright ' + str( year ) + ' Andrea Picco.  |')
		print('|   Url: www.apicco.github.io/trajectory_alignment/   |')
		print('|-----------------------------------------------------|')
	
	elif not printit :

		return version

	else :

		raise AttributeError('Please, if you want to print the header (printit = True) or if you want to return the verion number only (printit = False).')


def load_directory(path , pattern = '.txt' , sep = None , comment_char = '#' , dt = None , t_unit = '' , coord_unit = '' , intensity_normalisation = 'None' , fill_trajectory = True , **attrs ):

	"""
	load_directory(path , pattern = '.txt' , sep = None , comment_char = '#' , dt = None , t_unit = '' , coord_unit = '' , intensity_normalisation = 'None' , **attrs ):
	loads all the trajectories listed in 'path', which have the same 'pattern'.
	columns are separated by 'sep' (default is None: a indefinite number of 
	white spaces). Comments in the trajectory start with 'comment_char'.
	
	intensity_normalisation can be: 'None' (no normalisation, default), 'Integral' (normalise over the integral of the fluorescence intensity), 
	or 'Absolute' (normalise the fluorescence intensity values between 0 and 1)"

	**attrs is used to assign columns to the trajectory attributes and to 
	add annotations. 
	If the time interval is added (and 't' is not called in the **attrs) 
	then the time column 't' is added, and the 't_unit' can be set.
	If 'coord' is called then the unit must be added.
	fill_trajectory = True (default) fills missing frames with nan
	"""

	if ('coord' in attrs.keys()) & (len(coord_unit) == 0): 
		raise AttributeError('Please, specify the coordinate unit \'coord_unit\'')
	if ('t' in attrs.keys()) & (len(t_unit) == 0): 
		raise AttributeError('Please, specify the time unit \'t_unit\'')
	if (dt != None) & (len(t_unit) == 0): 
		raise AttributeError('Please, specify the time unit \'t_unit\'')
	if (dt != None) & ('t' in attrs.keys()):
		raise AttributeError('Time is already loaded by the trajectories, you cannot also compute it from frames. Please, either remove the dt option or do not load the \'t\' column from the trajectories')

	trajectories = [] #the list of trajectories
	if ( pattern[ len( pattern ) - 1 ] == '$' ) : 
		files = [ f for f in os.listdir(path) if f.endswith( pattern[ : - 1 ] ) ] #list all the files in path that have pattern
	else : 
		files = [ f for f in os.listdir(path) if pattern in f] #list all the files in path that have pattern

	print( 'Loading of trajectory files' )
	for file in files:

		trajectory = Traj(experiment = path, path = os.getcwd()+'/'+path)
		trajectory.load(path+'/'+file,sep = sep, comment_char = comment_char, **attrs)
		if (dt != None):
			trajectory.time(dt,t_unit)
		if ('coord' in attrs.keys()):

			trajectory.annotations('coord_unit',coord_unit)

		if intensity_normalisation == 'Integral' :
			
			trajectory.scale_f()

		elif intensity_normalisation == 'Absolute' :
		
			trajectory.norm_f()

		elif intensity_normalisation != 'None' :

			raise AttributeError( "load_directory: Please, choose a value for the variable intensity_normalisation between 'None' (no normalisation, default), 'Integral' (normalise over the integral of the fluorescence intensity), or 'Absolute' (normalise the fluorescence intensity values between 0 and 1)" )

		trajectory.annotations( 'intensity_normalisation' , intensity_normalisation )
		if fill_trajectory : trajectory.fill()
		trajectories.append(trajectory)
	
	print( "\n >> load_directory: The 'intensity_normalisation' applied to the trajectories is '" + intensity_normalisation + "' <<\n" )

	print( 'Loading of trajectory files ended' )
	return trajectories 

def MSD(input_t1 , input_t2):

	"""
	MSD(t1,t2): finds the rototranslation the  minimises the mean square displacement between the trajectories t1 and t2 and returns the rototranslation of t2.
	Adapted from Horn, 1987, to the 2D case with means weighted on the product of the fluorescence intensities.
	"""

	msdt1 = cp.deepcopy(input_t1)
	msdt2 = cp.deepcopy(input_t2)


	if (len(msdt1.f()) == 0) | (len(msdt2.f()) == 0): 
		raise AttributeError('MSD(msdt1,msdt2) requires that trajectories msdt1 and msdt2 have values for the fluorescence intensity')

	#the following code follow Horn's (1987) nomenclature. msdt1 is what is 
	#called in the paper as 'right coordinates'. 
	#msdt2 is what is called as 'left coordinates'
	
	with wr.catch_warnings():
		# if both f() are 0 or if their product is 0,  a warning about invalid true divide is output. Here we suppress such warnings.
		wr.simplefilter("ignore", category=RuntimeWarning)
		w = msdt1.f() * msdt2.f() / np.nansum( msdt1.f() * msdt2.f() )
	#computed the center of mass, weigthed on the fluorescence intensity product
	rc = np.array([ np.nansum( w * msdt1.coord()[0] ), np.nansum( w * msdt1.coord()[1] )])
	lc = np.array([ np.nansum( w * msdt2.coord()[0] ), np.nansum( w * msdt2.coord()[1] )])

	#translate the trajecotries to their weigthed center of mass
	msdt1.translate( -1 * rc )
	msdt2.translate( -1 * lc )

	Sxx = np.nansum( w * msdt2.coord()[0] * msdt1.coord()[0] )
	Sxy = np.nansum( w * msdt2.coord()[0] * msdt1.coord()[1] )
	Syx = np.nansum( w * msdt2.coord()[1] * msdt1.coord()[0] )
	Syy = np.nansum( w * msdt2.coord()[1] * msdt1.coord()[1] )

	A = ( Syx - Sxy )
	B = ( Sxx + Syy )

	if ( ( A == 0 ) & ( B == 0 ) ):
		theta = np.nan
	else : 
		theta = np.arctan2( - A , B )
#		M = - A / B
#
#		theta1 = np.arctan( M )
#		theta2 = theta1 + np.pi
#		
#		if B * np.cos(theta1) >= A * np.sin(theta1) :
#			theta = theta1
#		else:
#			theta = theta2

	msdt2.rotate( theta )

	#the 'score' is the mean square displacement weighted on the cross correlation of the fluorescence intensities
	score = np.nansum( w * ( msdt1.coord()[0] - msdt2.coord()[0] )**2 + w * ( msdt1.coord()[1] - msdt2.coord()[1] )**2 )
	
	#when the min_w is small and the software is sampling the start or end of trajectories,
	#which often have a large number of nan, then theta can become nan as M is the 
	#fraction of 0.0/0.0. If that happens then the score is set to infinite.
	if theta != theta: 
		score = np.inf

	return({ 
		'angle' : theta,
		'rc' : rc,
		'lc' : lc,
		'score' : score
		})

def nanMAD( x , axis = None , k = 1.4826):
	MAD = np.nanmedian( np.absolute( x - np.nanmedian( x , axis ) ) , axis )
	return( k * MAD )
	
def unified_start( t ) :

	try :

		return float( t.annotations()[ 'mean_starts' ] ) - 1.96 * float( t.annotations()[ 'std_starts' ] ) / np.sqrt( float( t.annotations()[ 'n_starts' ] ) )
		
	except :
		
		print( 'Error: one or more of the annotations mean_starts, std_starts, and n_starts is/are missiong' )

def unified_end( t ) :

	try :

		return float( t.annotations()[ 'mean_ends' ] ) + 1.96 * float( t.annotations()[ 'std_ends' ] ) / np.sqrt( float( t.annotations()[ 'n_ends' ] ) )
		
	except :
		
		print( 'Error: one or more of the annotations mean_ends, std_ends, and n_ends is/are missiong' )


def compute_average_start_and_end( trajectories_time_span , aligned_trajectories , max_frame ) :
	
	l = len( trajectories_time_span[ 'old_start' ] )
	# Define the START and the END of the average trajectory:
	#compute the start and the end of the average trajectory using the 
	#start and end of the aligned trajectories whose frame numbers are
	# greater than 0 (i.e. that appeared after the movie recording was 
	#started)
	traj_starts_to_average = [trajectories_time_span[ 'new_start' ][ j ] for j in range(l)\
			if trajectories_time_span[ 'old_start' ][ j ] > 0]
	if len( traj_starts_to_average ) > 0 : 
		mean_starts = np.mean( traj_starts_to_average )
		n_starts =  len( traj_starts_to_average )
		std_starts = np.std( traj_starts_to_average ) 
		
	else  :
		#it can be that all trajectories start with 0 (old_start), which means they 
		#started before the movie begun. If so the mean start is set as the latest 
		#time between the two trajectories.
		mean_starts = max(trajectories_time_span[ 'new_start' ])
		n_starts = np.nan
		std_starts = np.nan
		print( 'Warning: all trajectory starts were trunkated' )
	
	#same as for nan mean_starts. However, for the selected mean_ends is the smallest
	traj_ends_to_average = [trajectories_time_span[ 'new_end' ][ j ] for j in range(l)\
			if trajectories_time_span[ 'old_end' ][ j ] < ( max_frame - 3 ) * float(aligned_trajectories[ 0 ].annotations()[ 'delta_t' ])]
	if len( traj_ends_to_average ) > 0 : 
		mean_ends = np.mean( traj_ends_to_average )
		n_ends = len( traj_ends_to_average )
		std_ends = np.std( traj_ends_to_average ) 
	else  : 
		mean_ends = min(trajectories_time_span[ 'new_end' ])
		n_ends = np.nan
		std_ends = np.nan
		print( 'Warning: all trajectory ends were trunkated' )

	return( mean_starts , std_starts , n_starts , mean_ends , std_ends , n_ends )

#-------------------------------------END-OF-DEFINITIONS-in-compute_average_start_and_end-----------------------------------
def trajectory_average( aligned_trajectories_to_average , r , median , fimax ) :	

	#define the trajectory where the average will be stored
	t = Traj()

	#inherit the annotations from the reference trajectory
	for a in aligned_trajectories_to_average[ r ].annotations().keys():

		if a == 'file':
			t.annotations( 'reference_file' , aligned_trajectories_to_average[ r ].annotations()[ a ])
		else :
			t.annotations( a , aligned_trajectories_to_average[ r ].annotations()[ a ]) 

	if fimax :
		t.annotations( 'fimax' , 'TRUE' )

	#group all the attributes of the aligned trajectories...
	attributes = [ a for a in aligned_trajectories_to_average[ r ].attributes() if a not in ('t','frames')] 
	#create an empy dictionary where all the attributes that will be then averaged are stored
	attributes_to_be_averaged = {}
	for a in attributes:
		attributes_to_be_averaged[a] = []

	#merge all the trajectory attributes into the dictionary attributes_to_be_averaged.
	for j in range( len( aligned_trajectories_to_average ) ):
		
		for a in attributes:

			attributes_to_be_averaged[a].append(getattr(aligned_trajectories_to_average[ j ],'_'+a))

	#all the aligned trajectories are set to start at the same  mean_start and finish at mean_end computed from
	#trajectories_time_span in compute_average().Hence, the time interval is the same
	t.input_values( 't' , aligned_trajectories_to_average[ r ].t()) 
		
	#average the attributes of the trajectories and assign 
	#them to the average trajectory [ r ]
	with wr.catch_warnings():
		
		# if a line is made only of nan that are averaged, then a waring is outputed. Here we suppress such warnings.
		wr.simplefilter("ignore", category=RuntimeWarning)
	
		for a in attributes: 

			if a[ len( a ) - 4 : len( a ) ] == '_err' :
				
				raise AttributeError('The trajectories to be averaged have already an non empty error element, suggeting that they are already the result of an average. These error currently are not propagated. Check that your trajectories are correct')
			
			#if _a_err is in the trajectories slots, it means that the attribute a is not 
			#an error attribute (i.e. an attribute ending by _err; in fact if a would end by '_err'
			#then _a_err would have twice the appendix _err (i.e. _err_err) and would have 
			#no equivalent in the trajectory __slots__. If _a_err is then in the trajectory
			#slots, then both the mean and the sem can be computed. There is no sem without mean.

			if '_' + a + '_err' in t.__slots__:
				
				if median :

					t.input_values( a ,
							np.nanmedian( attributes_to_be_averaged[ a ], axis = 0 )
						)

				else :

					t.input_values( a ,
							np.nanmean( attributes_to_be_averaged[ a ], axis = 0 )
						)

				#if there is no n defined yet, it computes it
				if not t.n().any() : 
					#compute the number of not-nan data points by dividing
					#the nansum by the nanmean. The operation is performed
					#on the last attribute in the loop that can either have 
					#two dimensions (as 'coord') or one. In case of two dims
					#only one is used to compute '_n'.
					
					x = np.nanmean( attributes_to_be_averaged[ a ] , axis = 0 )
					y = np.nansum( attributes_to_be_averaged[a] , axis = 0 )
			
					if len(x) == 2:
						t.input_values( 'n' , y[0]/x[0] )
					else:
						with wr.catch_warnings():
							# if both y[i] and x[i] are 0, then a waring is outputed. Here we suppress such warnings.
							wr.simplefilter("ignore", category=RuntimeWarning)
							t.input_values( 'n' , y/x )

				#compute the errors as standard errors of the mean/median
				try :
					if median :

						t.input_values( a + '_err' ,
							nanMAD( attributes_to_be_averaged[ a ], axis = 0 ) / np.sqrt( t.n() )
							)

					else :

						t.input_values( a + '_err' ,
							np.nanstd( attributes_to_be_averaged[ a ], axis = 0 ) / np.sqrt( t.n() )
							)

				except :

					raise AttributeError( 'The attribute ' + a + ' cannot have its error assigned' )

			else :

				raise AttributeError( 'The attribute ' + a + ' is not recongnised as an attribute' )

	return( t )
#-------------------------------------END-OF-DEFINITION-of-trajectory_average-----------------------------------
def lie_down( t ):
	translation_vector = ( - np.nanmedian( t.coord()[0] ) ,- np.nanmedian( t.coord()[1] ) )
	t.translate( translation_vector )
	
	I_xx = np.nansum( t.f() * t.coord()[1] ** 2 )
	I_yy = np.nansum( t.f() * t.coord()[0] ** 2 )
	I_xy = np.nansum( t.f() * t.coord()[0] * t.coord()[1] )
	
	theta = np.arctan2( 2 * I_xy , I_xx - I_yy ) / 2

	I_x = I_xx + I_xy * np.tan( theta ) 
	I_y = I_yy - I_xy * np.tan( theta )

	if I_x > I_y : theta = theta - np.pi/2
	t.rotate( theta )

	with wr.catch_warnings():
		# if a coord has nan then a waring is outputed when nan > 0 or nan < 0 is asked. Here we suppress such warnings.
		wr.simplefilter("ignore", category=RuntimeWarning)
		A = np.nanmedian( t.coord()[0][ t.coord()[0] > 0 ] ** 2 )
		B = np.nanmedian( t.coord()[0][ t.coord()[0] < 0 ] ** 2 )

		if B > A : 
			t.rotate( np.pi )
			theta = theta + np.pi #to ouptput the angle
	
	model = linear_model.LinearRegression()	
	model_RANSACR = linear_model.RANSACRegressor( model , random_state = 42 )
	
	l = len( t )
	
	j = 0
	for i in range( l ) :
		
		if ( not np.isnan( t.coord()[ 0 ][ i ] ) ) & ( not np.isnan( t.coord()[ 1 ][ i ] ) ) :
			if j == 0 :
				X = np.array( [[ t.coord()[ 0 ][ i ] ]] )
				y = np.array( [ t.coord()[ 1 ][ i ] ] )
			else :
				X = np.insert( X , 0 , t.coord()[ 0 ][ i ] , axis = 0 )
				y = np.insert( y , 0 , t.coord()[ 1 ][ i ] , axis = 0 )
			j += 1

	with wr.catch_warnings():
		# also a bug warning occurs from linear models, RANSACR.
		wr.simplefilter("ignore", category=RuntimeWarning)
		model_RANSACR.fit( X , y )
	
	t.rotate( - np.arctan( model_RANSACR.estimator_.coef_[0] ) )

	return( { 'translation' :  translation_vector , 'angle' : theta - np.arctan( model_RANSACR.estimator_.coef_[0] ) } )


def average_trajectories( trajectory_list , output_file = 'average' , median = False , unify_start_end = False , max_frame=[] , fimax = False , fimax_filter = [ -3/35 , 12/35 , 17/35 , 12/35 , -3/35 ] ):

	"""
	average_trajectories( trajectory_list , max_frame = 500 , output_file = 'average' , median = False ): align all the 
	trajectories in the list together, and average them. 'max_frame' is the max number of frames the movies from which 
	the trajctories have been extracted has. It is used to check whether some trajectories are trunkated at the end.
	'output_file' is the name of the output. average_trajectories outputs a txt file with the average trajectroy and 
	a directory with all the raw trajectories that have been used to compute the average aligned together in space and time.
	median is an option to compute the median instead of the average of the aligned trajectories. It is useful in case 
	of noisy datasets.
	"""

	if len(trajectory_list) == 0 : 

		raise IndexError('There are not tajectories in the list; check that the trajectories were loaded correctly') 

	if not max_frame :

		raise TypeError('You need to specify the max_frame, which is the frame number in your movies')

	def R(alpha):
		"""
		R(alpha): returns the rotation matrix:
		
			/	cos(alpha)	-sin(alpha)		\
		R =	|								| 
			\	sin(alpha)	cos(alpha)		/
		"""
		return(np.matrix([[np.cos(alpha),-np.sin(alpha)], [np.sin(alpha),np.cos(alpha)]]))

	def float_range(x,y,step):
		while x < y:
			yield x
			x += step
	def compute_score(average_t,t_list):

		s = []
	
		for t in t_list:

			w = average_t.f() * t.f() / np.nansum( average_t.f() * t.f() )
			s.append(np.nansum(w * np.nansum((average_t.coord() - t.coord())**2)))
		
		return(s)

	def triplicate_trajectory( t ):
		#triplicate t adding itself at its beginning and at its end
	
		output = cp.deepcopy( t )
	
		#anticipate the start of trajectory by the trajectory duration and a time interval (you need one time interval
		#between the beginning of the real trajectory and the last point of the "anticipated" bit.
		#as len(t) is the number of frames + 1, then len(t) *dt is the duration of the trajectory + one dt interval, which is needed to
		#separate the duplicate trajectory form the original trajectory
		output.start( t.start() - ( len( t ) * float(t.annotations()['delta_t']) ) )
		#delay the end of trajectory in the same way 
		output.end( t.end() + ( len( t ) * float(t.annotations()['delta_t']) ) )
		
		#coord
		output.coord()[:,0:len(t)] = t.coord()
		output.coord()[:,( len(output) - len(t) ):len(output)] = t.coord()
		#f	
		output.f()[0:len(t)] = t.f()
		output.f()[( len(output) - len(t) ):len(output)] = t.f()

		return(output)
	def meanangle(angle_estimates):
		
		#angles can be identical +- n * pi. Hence, averaging 
		#absolute values for the mean angle is wrong. Example:
		# np.mean( ( 0, 2 * np.pi ) ) 
		#is not 0 but np.pi. However, cos and sin are invariant 
		#and the mean angle can be computed back from the 
		#mean cos and mean sin.
		#angle_estimates is an array (matrix) where the i-th
		#element has estimates of the angles that rotate the 
		#i-th trajectory.
		
		mean_cos = np.mean(np.cos(angle_estimates),axis=1)
		mean_sin = np.mean(np.sin(angle_estimates),axis=1)
	
		mean_angle = np.arctan2( mean_sin , mean_cos )
		
		return( mean_angle )

	def refine_alignment( t1 , t2 , lag , alignments , WeightTrajOverlap = False ):

		t1_frames = t1.frames()
		t2_frames = t2.frames() + lag

		sel_t1 = [ i for i in range( len(t1_frames) ) if t1_frames[i] in t2_frames ]
		sel_t2 = [ i for i in range( len(t2_frames) ) if t2_frames[i] in t1_frames ]
	
		if ( len(sel_t1) > 0 ) & ( len(sel_t2) > 0 ) :

			alignments.append(
					MSD( t1.extract( sel_t1 ) , t2.extract( sel_t2 ) )
					)
			alignments[ len( alignments ) - 1][ 'lag' ] = lag
	
			if WeightTrajOverlap :
				#the scores are weighted with the number of datapoints of the two trajectoreis that
				#trajectories (for example, trajectories that overlap with two data points only).
				
				alignments[ len( alignments ) - 1 ][ 'score' ] =\
						alignments[ len( alignments ) - 1 ][ 'score' ] / np.sqrt( len( sel_t1 ) )
		return()

		#-------------------------------------END-OF-DEFINITIONS-in-average_trajectories-----------------------------------

	def compute_transformations( t1 , t1_index , trajectory_list , fimax , fimax_filter ) :

			selected_alignments = []
			
			if ( fimax ) :
				
				t1 = t1.fimax( fimax_filter )
	
			for traj2 in [ traj2 for traj2 in trajectory_list ]:
				
				if trajectory_list.index(traj2) >= t1_index :
					
					#list of trajectories called in the second loop; as the transformation matrices are 
					#symmetric, transformations are computed only in the upper diagonal. 
					selected_alignments.append(
						{
							'angle' : 0,
							'rc' : np.array([0,0]),
							'lc' : np.array([0,0]),
							'lag' : 0,
							'lag_unit' : 'frames',
							'score' : np.NaN
							})
	
				else :

					t2 = cp.deepcopy( traj2 )
					#t2.norm_f()

					print( 'ref. traj.:\t' + t1.annotations()['file'] )
					print( 'aligned traj.:\t' + t2.annotations()['file'] )

					alignments = []
	
					if ( fimax ) : 
	
						t2 = t2.fimax( fimax_filter )
	
					#triplicate the longest trajectory by adding itself at its beginning and at its end
					if ( len(t1)  >= len(t2) ) :
						x = triplicate_trajectory(t1)
						y = t2
					else :
						x = triplicate_trajectory(t2)
						y = t1
	
					convolution_steps = len(x) - len(y) 
					for i in range( 0 , convolution_steps ) :
						#by triplicating the longest trajectory we can test all possible alignments in
						#space and time starting with the entire trajectories x and y.
						lag = int( x.frames( 0 ) - y.frames( 0 ) + i )
						x_frames = x.frames()
						y_frames = y.frames() + lag
						
						#select the frames that are overlapping 
						sel_frames = [ i for i in range( len( x_frames) ) if x_frames[ i ] in y_frames ]
			
						#which trajectory was triplicated decides the sign of the lag
						if ( len(t1)  >= len(t2) ) :
	
							alignments.append(
									MSD( x.extract( sel_frames ) , y )
									)
							alignments[ len(alignments)-1 ][ 'lag' ] = lag
						
						else :
					
							alignments.append(
									MSD( y , x.extract( sel_frames ) )
									)
							alignments[ len(alignments)-1 ][ 'lag' ] = - lag
					
					s = [ a['score'] for a in alignments ]
					lags = [ a['lag'] for a in alignments ]
					sel_alignments = [ i for i in range(len(s)) if s[i] == min(s) ]
				
					#check which of the selected alignments best fit the trajectory t1 
					#and not just its triplicate. Importantly, also recompute the alignment
					#without repetitions of the trajectory, which alter the alignment output
					refined_alignments_1 = []
					t1_frames = t1.frames()
					for sa in sel_alignments: 
	
						refine_alignment( t1 , t2 , lags[ sa ] , refined_alignments_1 , WeightTrajOverlap = True ) 
					
					refined_s_1 = [  a['score'] for a in refined_alignments_1 ]
					
					refined_alignments_2 =[]
					lag = refined_alignments_1[ refined_s_1.index( min( refined_s_1 ) ) ][ 'lag' ]
	
					#define a span, which is not too small, nor too big compared to the trajectory length
					refine_span = int( min( len( t1 ) , len( t2 ) ) / 10 )
					for refined_lag in range( lag - refine_span , lag + refine_span + 1 ):
	
						refine_alignment( t1 , t2 , refined_lag , refined_alignments_2 , WeightTrajOverlap = False )
	
					refined_s_2 = [  a['score'] for a in refined_alignments_2 ]
				
					selected_alignments.append( refined_alignments_2[ refined_s_2.index( min( refined_s_2 ) ) ] )
	
			if ( fimax ) :

				print('\nfimax = True; Transformations were computed using only the trajectory information up to the max in fluorescence intensity.')
	
			print('________________')
			
			return( selected_alignments )

	#-------------------------------------END-OF-DEFINITIONS-in-compute_transformations-----------------------------------

	def compute_average( trajectory_list , tranformations , median , fimax , max_frame , unify_start_end ) :
		
		aligned_trajectories = [] #contains all the alignments in respect to each trajectory
		average_trajectory = [] #contains all averages in respect to each trajectory
		alignment_precision = [] #contains the alignment precision, measured as a score of the alignment
	
		#As each trajectory is aligned to a reference trajectory or 
		#acts as a reference the rc and lc vectors are obtained 
		#from rcs and its transpose (i.e. the aligning trajectory
		#becomes the aligned trajectory).
		rcs = transformations['rcs'] + np.transpose(transformations['lcs'],axes=(1,0,2))

		l = len(transformations['angles'])
		#reference trajectories are indexed with r
		for r in range( l ) :
		
			#define a dictionary used to store the starts and ends of the aligned
			#trajectories to compute the start of the average trajectory
			trajectories_time_span = \
					{ 'old_start' : [], 'new_start' : [], 'old_end' : [], 'new_end' : []}
			
			#compute the transformation of the trajectories 
			#in respect to the r-th trajectory
			#--angles--
			angles_in_respect_of_r = transformations['angles'][ r , ] - transformations['angles']  
			m_angles = meanangle(angles_in_respect_of_r)
			#--lags--
			#lags_in_respect_of_r = transformations['lags'] - transformations['lags'][ r ,]
			lags_in_respect_of_r = transformations['lags'][ r ,] - transformations['lags']
			m_lags = [ int(round(l)) for l in np.mean(lags_in_respect_of_r,axis=1)]
			#--translations--
			r_cm = np.mean([rcs[ r , j ] for j in range(l) if j != r ] , axis = 0 )

			#make a copy of the trajectory_list, whose trajectories need to be aligned
			aligned_trajectories.append( cp.deepcopy( trajectory_list ) )
	
			##################################################	
			#align the trajectoris together in space and time
			##################################################	
			for j in range(l):
			
				trajectories_time_span[ 'old_start' ].append(aligned_trajectories[ r ][ j ].start())
				trajectories_time_span[ 'old_end' ].append(aligned_trajectories[ r ][ j ].end())
				
				#compute the center of mass of the full trajectory
	
				l_cm = np.mean([rcs[ j , r ] for r in range(l) if r != j ] , axis=0 )
		
				# the following is equivalent to
				#
				# R( m_angles ) @ aligned_trajectories + T
				#
				# where R would be the rotation matrix computed from m_angles
				# and T is the translation computed as
				#
				# r_cm - R( m_angles ) @ l_cm
				#
				# see Horn 1987 for details.
				aligned_trajectories[ r ][ j ].translate( - l_cm )
				aligned_trajectories[ r ][ j ].rotate( m_angles[ j ] )
	
				aligned_trajectories[ r ][ j ].translate( r_cm )
				aligned_trajectories[ r ][ j ].lag( m_lags[ j ] )

				aligned_trajectories[ r ][ j ].annotations()[ 'l_cm' ] = tuple( l_cm )
				aligned_trajectories[ r ][ j ].annotations()[ 'r_cm' ] = tuple( r_cm )
				aligned_trajectories[ r ][ j ].annotations()[ 'm_angle' ] = m_angles[ j ]
				aligned_trajectories[ r ][ j ].annotations()[ 'm_lag' ] = m_lags[ j ]
				
				trajectories_time_span[ 'new_start' ].append(aligned_trajectories[ r ][ j ].start())
				trajectories_time_span[ 'new_end' ].append(aligned_trajectories[ r ][ j ].end())

			mean_start , std_start , n_start , mean_end , std_end , n_end = compute_average_start_and_end( trajectories_time_span , aligned_trajectories[ r ] , max_frame )


			if unify_start_end :

				#uniform start and end of aligned trajectories to mean_start and mean_end
				for j in range(l):

					#record the standard deviation of the average start and end, so that the user
					#knows how the start and end timepoints of the raw trajectories are distributed,
					#once alingned.
					aligned_trajectories[ r ][ j ].annotations( 'mean_starts' , str( mean_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'std_starts' , str( std_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'n_starts' , str( n_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'mean_ends' , str( mean_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'std_ends' , str( std_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'n_ends' , str( n_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'unify_start_end' , str( unify_start_end ) )
		
					#if the unify_start_end is choosen, the average trajectory is started (and ended) from the average 
					#start (and end) of the trajectory minus (and plus) the 95% CI. This addition (or subtraction) as
					#been choosen to counter the intrinsic underestimate of the trajectories lifetimes.
					aligned_trajectories[ r ][ j ].start( unify_start( aligned_trajectories[ r ][ i ] ) )
					aligned_trajectories[ r ][ j ].end( unify_end( aligned_trajectories[ r ][ i ] ) )

			else :

				for j in range(l):
	
					#record the standard deviation of the average start and end, so that the user
					#knows how the start and end timepoints of the raw trajectories are distributed,
					#once alingned.
					aligned_trajectories[ r ][ j ].annotations( 'mean_starts' , str( mean_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'std_starts' , str( std_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'n_starts' , str( n_start ) )
					aligned_trajectories[ r ][ j ].annotations( 'mean_ends' , str( mean_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'std_ends' , str( std_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'n_ends' , str( n_end ) )
					aligned_trajectories[ r ][ j ].annotations( 'unify_start_end' , str( unify_start_end ) )
		
					aligned_trajectories[ r ][ j ].start( min( trajectories_time_span[ 'new_start' ] ) )
					aligned_trajectories[ r ][ j ].end( max( trajectories_time_span[ 'new_end' ] ) )
	
			########################################################################	
			#compute the average of the trajectories aligned to the r-th trajectory
			#define the average trajectory and its time attribute
			########################################################################	
		
			ta = trajectory_average( aligned_trajectories[ r ] , r , median , fimax ) 
			
#TO DEL			if not unify_start_end :
#TO DEL
#TO DEL				ta.annotations( 'unified_start' , mean_start - 1.96 * std_start / np.sqrt( n_start ) )
#TO DEL				ta.annotations( 'unified_end' , mean_end + 1.96 * std_end / np.sqrt( n_end ) )
			
			average_trajectory.append( ta )

			#store the transformations of the trajectories in respect of the trajectory r.
			if r == 0:
				all_m_angles = np.array([ m_angles ])
				all_m_lags = np.array([ m_lags ])
			else :
				all_m_angles = np.vstack([ all_m_angles , m_angles ])
				all_m_lags = np.vstack([ all_m_lags , m_lags ])
	
			# make a copy of the average trajectory ta, and unify its start and end 
			# to compute a mean precision that reflects the invagination dynamice
			# and not how well noisy and/or excessively long trajectories might
			# align.
			ta_tmp = cp.deepcopy( ta )
			ta_tmp.start( unified_start( ta_tmp ) )
			ta_tmp.end( unified_end( ta_tmp ) )

			mean_precision =  np.sqrt(
					np.nanmean( 
						ta_tmp.coord_err()[ 0 ] ** 2 + ta_tmp.coord_err()[ 1 ] ** 2 
						)
					)
			alignment_precision.append(mean_precision)
		
		print('ALIGNMENT PRECISIONS.\nMIN is the alignment\nselected for the average\n----------------------')
		for a in alignment_precision :
			
			if a == min( alignment_precision ) :

				print( 'MIN>>\t' + str( a ) )
			
			elif a == max( alignment_precision ) :
		
				print( 'MAX>>\t' + str( a ) )
			
			else :
		
				print( '\t' + str( a ) )

		print('----------------------')
		print( 'MEAN:\t' + str( np.nanmean( alignment_precision ) ) )

		return( aligned_trajectories , average_trajectory , alignment_precision )
	
	#-------------------------------------END-OF-DEFINITIONS-in-compute_average-----------------------------------

	header() 

	print( '\nunify_start_end = ' + str( unify_start_end ) )
	
	#define the list where transformations are stored
	transformations = {
			'angles' : np.array( [] ),
			'rcs' : np.array( [ np.array( [] ) , np.array( [] ) ] ),#note that the matric rcs is the transpose of the lcs
			'lcs' : np.array( [ np.array( [] ) , np.array( [] ) ] ),
			'lags' : np.array( [] ),
			'lag_units' : np.array( [] )
			}

	for traj1 in trajectory_list:

		#t1 is the reference trajectory to which all the other trajectories are alinged
		#The loop goes on all trajectories as all of them are eligible to be used as reference

		
		#the index need to be computed now, becuase if fi_max is true, then t1 will be replaced 
		#by the part of t1 trajectory that stops at the peak of fluorescence intensity. This 
		#new trajectory cannot be found animore in trajectory_list. Hence, we must compute the 
		#index before.
		
		t1 = cp.deepcopy( traj1 )
		t1_index = trajectory_list.index(traj1) 	
		#t1.norm_f()

		selected_alignments = compute_transformations( t1 , t1_index , trajectory_list , fimax , fimax_filter )

		#Create a matrix with all the transformations: angle, lag and center of masses. 
		#As a convention the element i,j in the matrix contains the elements for the
		#rototranslation and temporal shift to align the trajectori i to j, j being the
		#reference.
		if t1_index == 0:

			transformations['angles'] = np.array(
					[a['angle'] for a in selected_alignments]
					)
			transformations['rcs'] = np.array([
					np.array([a['rc'] for a in selected_alignments])
					])
			transformations['lcs'] = np.array([
					np.array([a['lc'] for a in selected_alignments])
					])
			transformations['lags'] = np.array(
					[a['lag'] for a in selected_alignments]
					)
		else:

			transformations['angles'] = np.vstack([
				transformations['angles'],
				np.array([
					np.array([a['angle'] for a in selected_alignments])
					])
				])
			transformations['rcs'] = np.vstack([
				transformations['rcs'],
					[[a['rc'] for a in selected_alignments]]
				])
			transformations['lcs'] = np.vstack([
				transformations['lcs'],
					[[a['lc'] for a in selected_alignments]]
				])
			transformations['lags'] = np.vstack([
				transformations['lags'],
				np.array(
					[a['lag'] for a in selected_alignments]
					)
				])

	transformations['angles'] = transformations['angles'] - np.transpose(transformations['angles'])
	transformations['lags'] = transformations['lags'] - np.transpose(transformations['lags'])
	
	l = len(transformations['angles'])
	for i in range( l ):
		transformations['lcs'][ i , i ] = [ 0 , 0 ]
	
	#compute the average transformation using each trajectory as possible reference
	aligned_trajectories , average_trajectory , alignment_precision = compute_average( trajectory_list , transformations , median , fimax , max_frame , unify_start_end )

	best_average = alignment_precision.index( min( alignment_precision ) ) 
	worst_average = alignment_precision.index( max( alignment_precision ) ) 

	if not unify_start_end :

		# compute the lie_down only on the part of the trajectory that represents
		# most of the average trajectories. That would be the part of average trajectory 
		# chosen if unify_start_end = True, i.e. the part of trajectory comprised between
		# the annotations unified_start and unified_end
		average_trajectory_tmp = cp.deepcopy( average_trajectory[ best_average ] )
		average_trajectory_tmp.start( unified_start( average_trajectory_tmp ) )
		average_trajectory_tmp.end( unified_end( average_trajectory_tmp ) )
#TO DEL		average_trajectory_tmp.start( float( average_trajectory_tmp.annotations()[ 'unified_start' ] ) )
#TO DEL		average_trajectory_tmp.end( float( average_trajectory_tmp.annotations()[ 'unified_end' ] ) )
		lie_down_transform = lie_down( average_trajectory_tmp )

	else :
	
		lie_down_transform = lie_down( average_trajectory[ best_average ] )

	# lie_down transformations applied to average_trajectory[ best_average ]
	average_trajectory[ best_average ].translate( lie_down_transform[ 'translation' ] )
	average_trajectory[ best_average ].rotate( lie_down_transform[ 'angle' ] )

	average_trajectory[ best_average ].annotations()[ 'trajalign_version' ] = header( printit = False )
	average_trajectory[ best_average ].save( output_file )

	#save the trajectories use to compute the average, lied down as the average trajectory
	for i in range(l):

		aligned_trajectories[ best_average ][ i ].translate( lie_down_transform[ 'translation' ] )
		aligned_trajectories[ best_average ][ i ].rotate( lie_down_transform[ 'angle' ] )
		aligned_trajectories[ best_average ][ i ].annotations()[ 'trajalign_version' ] = header( printit = False )
		aligned_trajectories[ best_average ][ i ].annotations()[ 'lie_down_angle' ] = lie_down_transform[ 'angle' ]
		aligned_trajectories[ best_average ][ i ].annotations()[ 'lie_down_translation' ] = tuple( lie_down_transform[ 'translation' ] )

		filename = "./" + output_file + "/" + aligned_trajectories[ best_average ][ i ].annotations()[ 'file' ]
		if i == 0 :
			directory = os.path.dirname( filename )
			if not os.path.exists( directory ) :
				os.makedirs( directory )
		aligned_trajectories[ best_average ][ i ].save( filename )
	
	with open( "./" + output_file + "/alignment_precision.txt" , 'w' ) as f :

		for ap in alignment_precision :	
			f.write( repr( ap ) + '\n' )
	
	f.close()

	return( average_trajectory[ best_average ] , average_trajectory[ worst_average ] , { 'best_score' : aligned_trajectories[ best_average ] , 'worst_score' : aligned_trajectories[ worst_average ] } )

