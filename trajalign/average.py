#module of functions required to load and average the trajectories together
import os 
from trajalign.traj import Traj
import copy as cp
import numpy as np
import warnings as wr

from sklearn import linear_model

#import matplotlib
#matplotlib.use('Agg')
#from matplotlib import pyplot as plt

def load_directory(path,pattern = '.txt',sep = None,comment_char = '#',dt=None,t_unit='',coord_unit='',**attrs):

	"""
	load_directory(path,pattern = '.txt',sep = None,comment_char = '#',dt=None,t_unit='',**attrs)
	loads all the trajectories listed in 'path', which have the same 'pattern'.
	columns are separated by 'sep' (default is None: a indefinite number of 
	white spaces). Comments in the trajectory start with 'comment_char'.
	**attrs is used to assign columns to the trajectory attributes and to 
	add annotations. 
	If the time interval is added (and 't' is not called in the **attrs) 
	then the time column 't' is added, and the 't_unit' can be set.
	If 'coord' is called then the unit must be added.
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

	for file in files:
		trajectory = Traj(experiment = path, path = os.getcwd()+'/'+path, file = file)
		trajectory.load(path+'/'+file,sep = sep, comment_char = comment_char, **attrs)
		if (dt != None):
			trajectory.time(dt,t_unit)
		if ('coord' in attrs.keys()):
			trajectory.annotations('coord_unit',coord_unit)
		trajectory.fill()
		trajectory.norm_f()
		trajectories.append(trajectory)
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

def average_trajectories( trajectory_list , max_frame=500 , output_file = 'average' , median = False ):

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
	
	def R(alpha):
		"""
		R(alpha): returns the rotation matrix:
		
			/	cos(alpha)	-sin(alpha) \
		R =	|							| 
			\	sin(alpha)	cos(alpha)	/
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
	def triplicate_trajectory(t):
		#triplicate t adding itself at its beginning and at its end
	
		output = cp.deepcopy(t)
		
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

		#X = np.nanmedian( t.coord()[0] ** 2 )  - np.nanmedian( t.coord()[0] ) **2 
		#Y = np.nanmedian( t.coord()[0] * t.coord()[1] ) - np.nanmedian( t.coord()[0] ) * np.nanmedian( t.coord()[1] )
		#c = np.nanmedian( t.coord()[1] ) - Y * np.nanmedian( t.coord()[0] ) / X
	
		#t.translate( ( 0 , c ) )
		#t.rotate( np.arctan2( Y , X ) )

		with wr.catch_warnings():
			# if a coord has nan then a waring is outputed when nan > 0 or nan < 0 is asked. Here we suppress such warnings.
			wr.simplefilter("ignore", category=RuntimeWarning)
			A = np.nanmedian( t.coord()[0][ t.coord()[0] > 0 ] ** 2 )
			B = np.nanmedian( t.coord()[0][ t.coord()[0] < 0 ] ** 2 )

			if B > A : 
				t.rotate( np.pi )
				theta = theta + np.pi #to ouptput the angle
		
		model = linear_model.LinearRegression()	
		model_RANSACR = linear_model.RANSACRegressor( model )
		
		l = len( t )
		
		for i in range( l ) :
			if ( not np.isnan( t.coord()[ 0 ][ i ] ) ) & ( not np.isnan( t.coord()[ 1 ][ i ] ) ) :
				if i == 0 :
					X = np.array( [[ t.coord()[ 0 ][ i ] ]] )
					y = np.array( [ t.coord()[ 1 ][ i ] ] )
				else :
					X = np.insert( X , 0 , t.coord()[ 0 ][ i ] , axis = 0 )
					y = np.insert( y , 0 , t.coord()[ 1 ][ i ] , axis = 0 )

		with wr.catch_warnings():
			# also a bug warning occurs from linear models, RANSACR.
			wr.simplefilter("ignore", category=RuntimeWarning)
			model_RANSACR.fit( X , y )
		
		t.rotate( - np.arctan( model_RANSACR.estimator_.coef_[0] ) )
		return( { 'translation' :  translation_vector , 'angle' : theta - np.arctan( model_RANSACR.estimator_.coef_[0] ) } )

	#NOTE: the R code computed the time aligment by the CC of the FI filtered with a running filter of length 5 if FIMAX = TRUE. The score of the alignments was also fitered.
	#If FIMAX = FALSE the score only was filtered.
	transformations = {
			'angles' : np.array([]),
			'rcs' : np.array([np.array([]),np.array([])]),#note that the matric rcs is the transpose of the lcs
			'lcs' : np.array([np.array([]),np.array([])]),
			'lags' : np.array([]),
			'lag_units' : np.array([])
			}
	for t1 in trajectory_list: 

		#t1 is the reference trajectory to which all the other trajectories are alinged
		#The loop goes on all trajectories as all of them are eligible to be used as reference

		selected_alignments = []

		for t2 in [ t2 for t2 in trajectory_list]:
			if trajectory_list.index(t2) >= trajectory_list.index(t1):
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
				print( 'ref. traj.:\t' + t1.annotations()['file'] )
				print( 'aligned traj.:\t' + t2.annotations()['file'] )
				alignments = []

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


#				plt.figure()
#				plt.subplot(211)
#				plt.plot( t1.frames() , t1.f() , '-' )
#				plt.plot( t2.frames() + lag , t2.f() , '-' )
#				plt.plot( t2.frames() + refined_alignments_2[ refined_s_2.index( min( refined_s_2 ) ) ]['lag'] , t2.f() , '-' )
#				plt.subplot(212)
#				plt.plot( lags , s )
#				plt.plot( [ a['lag'] for a in refined_alignments_2 ] , refined_s_2 )
#				plt.show()

				
		print('________________')

		#Create a matrix with all the transformations: angle, lag and center of masses. 
		#As a convention the element i,j in the matrix contains the elements for the
		#rototranslation and temporal shift to align the trajectori i to j, j being the
		#reference.
		if trajectory_list.index(t1) == 0:

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

	transformations['angles'] = transformations['angles']-np.transpose(transformations['angles'])
	transformations['lags'] = transformations['lags']-np.transpose(transformations['lags'])

	l = len(transformations['angles'])
	for i in range(l):
		transformations['lcs'][ i , i ] = [ 0 , 0 ]
	#As each trajectory is aligned to a reference trajectory or 
	#acts as a reference the rc and lc vectors are obtained 
	#from rcs and its transpose (i.e. the aligning trajectory
	#becomes the aligned trajectory).
	rcs = transformations['rcs'] + np.transpose(transformations['lcs'],axes=(1,0,2))

	#transformations
	print('--angles--')
	print(transformations['angles'])
	print('--lags--')
	print(transformations['lags'])
	print('full rcs')
	print(rcs)

	#compute the average transformation using each trajectory as possible reference
	aligned_trajectories = [] #contains all the alignments in respect to each trajectory
	average_trajectory = [] #contains all averages in respect to each trajectory
	alignment_precision = [] #contains the alignment precision, measured as a score of the alignment

	#reference trajectories are indexed with r
	for r in range(l) :
	
		#define a dictionary used to store the starts and ends of the aligned
		#trajectories to compute the start of the average trajectory
		trajectory_time_span = \
				{ 'old_start' : [], 'new_start' : [], 'old_end' : [], 'new_end' : []}\
		
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
		#aligned_trajectories.insert( r ,cp.deepcopy(trajectory_list))
		##################################################	
		#align the trajectoris together in space and time
		##################################################	
		for j in range(l):
		
			trajectory_time_span[ 'old_start' ].append(aligned_trajectories[ r ][ j ].start())
			trajectory_time_span[ 'old_end' ].append(aligned_trajectories[ r ][ j ].end())
			
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
			
			trajectory_time_span[ 'new_start' ].append(aligned_trajectories[ r ][ j ].start())
			trajectory_time_span[ 'new_end' ].append(aligned_trajectories[ r ][ j ].end())
			
#			if ( j != r ):
#				plt.plot(trajectory_list[ r ].coord()[ 0 ] ,trajectory_list[ r ].coord()[ 1 ] , 'r' )
#				plt.plot(aligned_trajectories[ r ][ j ].coord()[ 0 ],aligned_trajectories[ r ][ j ].coord()[ 1 ] , 'g' )
#				plt.plot(trajectory_list[ r ].t(),trajectory_list[ r ].f() , 'r' )
#				plt.plot(aligned_trajectories[ r ][ j ].t(),aligned_trajectories[ r ][ j ].f() , 'g' )
#				plt.show()

		########################################################################	
		#compute the average of the trajectories aligned to the r-th trajectory
		#define the average trajectory and its time attribute
		########################################################################	
		average_trajectory.append( Traj() )
	
		#inherit the annotations from the reference trajectory
		for a in aligned_trajectories[ r ][ r ].annotations().keys(): #IT WAS [ r ][ 0 ]
			if a == 'file':
				average_trajectory[ r ].annotations( 'reference_file' , aligned_trajectories[ r ][ r ].annotations()[ a ])
			else :
				average_trajectory[ r ].annotations( a , aligned_trajectories[ r ][ r ].annotations()[ a ]) #IT WAS [ r ][ 0 ]

		# Define the START and the END of the average trajectory:
		#compute the start and the end of the average trajectory using the 
		#start and end of the aligned trajectories whose frame numbers are
		# greater than 0 (i.e. that appeared after the movie recording was 
		#started)
		traj_starts_to_average = [trajectory_time_span[ 'new_start' ][ j ] for j in range(l)\
				if trajectory_time_span[ 'old_start' ][ j ] > 0]
		if len( traj_starts_to_average ) > 0 : 
			mean_start = np.mean( traj_starts_to_average )
		else  :
			#it can be that all trajectories start with 0 (old_start), which means they 
			#started before the movie begun. If so the mean start is set as the latest 
			#time between the two trajectories.
			mean_start = max(trajectory_time_span[ 'new_start' ])
			average_trajectory[ r ].annotations( 'Warning_start' , 'all trajectory starts were trunkated' )
		
		#same as for nan mean_start. However, for the selected mean_end is the smallest
		traj_ends_to_average = [trajectory_time_span[ 'new_end' ][ j ] for j in range(l)\
				if trajectory_time_span[ 'old_end' ][ j ] < ( max_frame - 3 ) * float(aligned_trajectories[ r ][ 0 ].annotations()[ 'delta_t' ])]
		if len( traj_ends_to_average ) > 0 : 
			mean_end = np.mean( traj_ends_to_average )
		else  : 
			mean_end = min(trajectory_time_span[ 'new_end' ])
			average_trajectory[ r ].annotations( 'Warning_end' , 'all trajectory ends were trunkated' )

		#group all the attributes of the aligned trajectories...
		attributes = [ a for a in aligned_trajectories[ r ][ r ].attributes() if a not in ('t','frames')] 
		#create an empy dictionary where all the attributes that will be then averaged are stored
		attributes_to_be_averaged = {}
		for a in attributes:
			attributes_to_be_averaged[a] = []

		#merge all the trajectory attributes into the dictionary attributes_to_be_averaged.
		for j in range(l):
			aligned_trajectories[ r ][ j ].start( mean_start )
			aligned_trajectories[ r ][ j ].end( mean_end )
			for a in attributes:
				attributes_to_be_averaged[a].append(getattr(aligned_trajectories[ r ][ j ],'_'+a))
	
		#all the aligned trajectories are set to start at mean_start and finish at mean_end
		#setattr(average_trajectory[ r ],'_t',aligned_trajectories[ r ][ r ].t()) #IT WAS [ r ][ j ]
		average_trajectory[ r ].input_values( 't' , aligned_trajectories[ r ][ r ].t()) 
			
		#average the attributes of the trajectories and assign 
		#them to the average trajectory [ r ]
		with wr.catch_warnings():
			
			# if a line is made only of nan that are averaged, then a waring is outputed. Here we suppress such warnings.
			wr.simplefilter("ignore", category=RuntimeWarning)
		
			for a in attributes: 

				#if _a_err is in the trajectories slots, it means that the attribute a is not 
				#an error attribute (i.e. an attribute ending by _err; in fact if a would end by '_err'
				#then _a_err would have twice the appendix _err (i.e. _err_err) and would have 
				#no equivalent in the trajectory __slots__. If _a_err is then in the trajectory
				#slots, then both the mean and the std can be computed. There is no std without mean.
	
				if a[len(a)-4:len(a)] == '_err' :
					
					raise AttributeError('The trajectories to be averaged have already an non empty error element, suggeting that they are already the result of an average. These error currently are not propagated. Check that your trajectories are correct')
				
				if '_' + a in average_trajectory[ r ].__slots__:

					if median :

						average_trajectory[ r ].input_values( a ,
								np.nanmedian( attributes_to_be_averaged[ a ], axis = 0 )
							)

					else :

						average_trajectory[ r ].input_values( a ,
								np.nanmean( attributes_to_be_averaged[ a ], axis = 0 )
							)

					try :
						if median :

							average_trajectory[ r ].input_values( a + '_err' ,
								nanMAD( attributes_to_be_averaged[ a ], axis = 0 )
								)

						else :

							average_trajectory[ r ].input_values( a + '_err' ,
								np.nanstd( attributes_to_be_averaged[ a ], axis = 0 )
								)

					except :

						raise AttributeError( 'The attribute ' + a + ' cannot have its error assigned' )

				else :

					raise AttributeError( 'The attribute ' + a + ' is not recongnised as an attribute' )

		#compute the number of not-nan data points by dividing
		#the nansum by the nanmean. The operation is performed
		#on the last attribute in the loop that can either have 
		#two dimensions (as 'coord') or one. In case of two dims
		#only one is used to compute '_n'.
		
		x = np.nanmean( attributes_to_be_averaged[ a ] , axis = 0 )
		y = np.nansum( attributes_to_be_averaged[a] , axis = 0 )

		if len(x) == 2:
			#setattr(average_trajectory[ r ],'_n',y[0]/x[0])
			average_trajectory[ r ].input_values( 'n' , y[0]/x[0] )
		else:
			with wr.catch_warnings():
				# if both y[i] and x[i] are 0, then a waring is outputed. Here we suppress such warnings.
				wr.simplefilter("ignore", category=RuntimeWarning)
				#setattr(average_trajectory[ r ],'_n',y/x)
				average_trajectory[ r ].input_values( 'n' , y/x )

		#store the transformations of the trajectories in respect of the trajectory r.
		if r == 0:
			all_m_angles = np.array([ m_angles ])
			all_m_lags = np.array([ m_lags ])
		else:
			all_m_angles = np.vstack([ all_m_angles , m_angles ])
			all_m_lags = np.vstack([ all_m_lags , m_lags ])
		
		mean_precision =  np.sqrt(
				np.nanmean( 
					average_trajectory[ r ].coord_err()[ 0 ] ** 2 + average_trajectory[ r ].coord_err()[ 1 ] ** 2 
					)
				)
		alignment_precision.append(mean_precision)

	best_average = alignment_precision.index( min( alignment_precision ) ) 
	worst_average = alignment_precision.index( max( alignment_precision ) ) 
	lie_down_transform = lie_down( average_trajectory[ best_average ] )
	average_trajectory[ best_average ].save( output_file )

	print('alignment_precision')
	print(alignment_precision)
	print('angles')
	print(all_m_angles)
	print('lags')
	print(all_m_lags)
	print(all_m_lags - all_m_lags[0])

	#save the trajectories use to compute the average, lied down as the average trajectory
	for i in range(l):
		aligned_trajectories[ best_average ][ i ].translate( lie_down_transform[ 'translation' ] )
		aligned_trajectories[ best_average ][ i ].rotate( lie_down_transform[ 'angle' ] )
		filename = "./" + output_file + "/" + aligned_trajectories[ best_average ][ i ].annotations( 'file' )
		if i == 0 :
			directory = os.path.dirname( filename )
			if not os.path.exists( directory ) :
				os.makedirs( directory )
		aligned_trajectories[ best_average ][ i ].save( filename )

	return( average_trajectory[ best_average ], average_trajectory[ worst_average] , aligned_trajectories[ best_average ] )
