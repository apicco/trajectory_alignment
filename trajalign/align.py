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
from trajalign.average import header , load_directory , nanMAD , MSD , compute_average_start_and_end , trajectory_average 
from scipy.interpolate import UnivariateSpline #need to install py35-scikit-learn
import numpy as np
import copy as cp

from matplotlib import pyplot as plt

def spline( t1 , t2 ) :

	"""
	interpolate t1 or t2 with a spline.
	"""

	#the interpolation function
	def interpolation( to_interpolate , delta_t , k = 3 ) :
	
		interpolated_traj = Traj( interpolated = 'True' )
		interpolated_traj.annotations( to_interpolate.annotations() )
		interpolated_traj.annotations()[ 'delta_t' ] = delta_t

		l = len( to_interpolate )

		if not l > k :
			
			#UnivariateSpline requires that m > k, where m is the number of points interpolated 
			#and k is the degree of smoothing spline. Default in UnivatiateSpline and in here is k = 3.
			k = l - 1

		#the new time intervals for the trajectory interpolated
		t = [ to_interpolate.start() ]
		while( t[ len(t) - 1 ] <= to_interpolate.end() ) :
			t.append( t[ len(t) - 1 ] + delta_t )
		
		interpolated_traj.input_values( 
				name = 't' , 
				x = t ,
				unit = to_interpolate.annotations( 't_unit' ) 
				)

		for attribute in to_interpolate.attributes() : 	
			
			if attribute in [ 'f' , 'mol' ] :

				s = UnivariateSpline( to_interpolate.t() , getattr( to_interpolate , '_'+attribute ) , k = k )
				interpolated_traj.input_values( 
						name = attribute , 
						x = s( interpolated_traj.t() )
						)

			if attribute == 'coord' :

				s_x = UnivariateSpline( to_interpolate.t() , to_interpolate.coord()[ 0 ] , k = k )
				s_y = UnivariateSpline( to_interpolate.t() , to_interpolate.coord()[ 1 ] , k = k )
				interpolated_traj.input_values( 
						name = 'coord' , 
						x = [ s_x( interpolated_traj.t() ) , s_y( interpolated_traj.t() ) ],
						)

		return( interpolated_traj )

	#the trajectory with the largest delta_t will be the one that will 
	#be splined. 

	if t1.annotations()[ 'delta_t' ] >= t2.annotations()[ 'delta_t' ] :

		delta_t = float(t2.annotations()[ 'delta_t' ])
	
	else :
		
		delta_t = float(t1.annotations()[ 'delta_t' ])

	not_nan = [ i for i in range( len( t1 ) ) if t1.f( i ) == t1.f( i ) ]
	t1_to_interpolate = t1.extract( not_nan )

	not_nan = [ i for i in range( len( t2 ) ) if t2.f( i ) == t2.f( i ) ]
	t2_to_interpolate = t2.extract( not_nan )

	return( 
			interpolation( t1_to_interpolate , delta_t ) ,
			interpolation( t2_to_interpolate , delta_t )
			)

def cc( input_t1 , input_t2 ):
	
	"""
	cc( input_t1 , input_t2 ) returns the time lag between the trajectory input_t1 and the trajectory input_t2,
	computed from the cross correlation of the fluorescence intensities of the two trajectories. 
	The trajectory input_t2 will be aligned in time to input_t1 by adding the output of cc to input_t2.t()
	"""

	t1 = cp.deepcopy( input_t1 )
	t2 = cp.deepcopy( input_t2 )

	if t1.annotations()[ 'delta_t' ] != t2.annotations()[ 'delta_t' ] :
		raise AttributeError('The two trajectories have different \'delta_t\' ') 
	else: 
		delta_t = t1.annotations()[ 'delta_t' ]
	
	#extend t1 to be as long as to include the equivalent
	#of t2 lifetime as NA before and after it:

	t1.start( t1.start() - t2.lifetime() )
	t1.end( t1.end() + t2.lifetime() )
	
	#align the two trajectories to the same start point
	lag0 = t1.start() - t2.start()
	t2.input_values( 't' , t2.t() + lag0 )

	output = []
	while t2.end() <= t1.end() :

		#because of rounding errors I cannot use:
		#f1 = [ t1.f( i ) for i in range( len( t1 ) ) if ( t1.t( i ) >= t2.start() ) & ( t1.t( i ) <= t2.end() ) ] 
		#but the following, where instead of greater than... I use >< delta_t / 2
		f1 = [ t1.f( i ) for i in range( len( t1 ) ) if ( t1.t( i ) - t2.start() > - delta_t / 2 ) & ( t1.t( i ) - t2.end() < delta_t / 2 ) ] 
		
		if len( f1 ) != len( t2 ) : raise IndexError( "There is a problem with the selection of t1 fluorescence intensities and t2 length in the cross-correlation function cc. The lengths do not match.")

		output.append( 
				sum( [ f1[ i ] * t2.f( i ) for i in range( len( t2 ) ) if ( f1[ i ] == f1[ i ] ) & ( t2.f( i ) == t2.f( i ) ) ] )
				)
		
		t2.lag( 1 )

	return( lag0 + output.index( max( output ) ) * t1.annotations()[ 'delta_t' ] )

def unify_start_and_end( t1 , t2 ):

	"""
	Uniform the start and the end of two trajectories that overlap
	in time, so that the overlapping time points can be used to compute the
	rotation and translation that aligns the two trajectories together.
	"""
	
	if t1.annotations()[ 'delta_t' ] != t2.annotations()[ 'delta_t' ] : 
		raise AttributeError('The trajectoires inputed in unify_start_and_end \
				have different delta_t')
	if t1.start() >= t2.end() : 
		raise AttributeError('The trajectory t1 inputed in unify_start_and_end \
				starts after the trajectory t2. The two trajectories must significantly overlap')
	if t2.start() >= t1.end() : 
		raise AttributeError('The trajectory t2 inputed in unify_start_and_end \
				starts after the trajectory t2. The two trajectories must significantly overlap')

	if t1.start() < t2.start() :
		t1.start( t2.start() )
	else :
		t2.start( t1.start() )
	if t1.end() < t2.end() :
		t2.end( t1.end() )
	else :
		t1.end( t2.end() )

	return()

def R( angle ) : 

	return( np.matrix( [[ np.cos( angle ) , - np.sin( angle ) ] , [ np.sin( angle ) , np.cos( angle ) ]] , dtype = 'float64' ) )

#-------------------------END-OF-DEFINITIONS--------------------------------

def align( path_target , path_reference , ch1 , ch2 , fimax1 = False , fimax2 = False , fimax_filter = [ -3/35 , 12/35 , 17/35 , 12/35 , -3/35 ] ):

	"""
	align( path_target , path_reference , ch1 , ch2 , ):
	aligns in space and in time the trajectories identified by path_target to path_reference,
	which is the reference trajectory. As a convention within the align function trajectories 
	labeled with 1 are the target trajectories that need to be ligned to the reference 
	trajectories, which are labelled with 2.The alignment uses the trajectories in ch1 
	and ch2, which have been acquired simultaneously and whose alignment has been 
	corrected for chormatic aberrations and imaging misalignments. 'ch1' refers to 
	the trajectories that need to be aligned to the average trajectory in 'path_target'. 
	'ch2' refers to 'path_reference'. Both the target and the reference trajectories can be
	aligned using only the trajectory information up to the peak of fluorescence intensity by 
	setting fimax1 and fimax2 to True, respectively. If fimax1 and/or fimax2 are true, then fimax_filer
	is used to compute where the peak of fluorescence intensity is. If no filter is desired, set 
	fimax_filer = [ 1 ].
	"""

	header() 

	#################################################################################################################
	#average trajectories are centered on their center of mass and must have been previously lied down 
	#(lie_down function in trajalign/average.py) so that they are orientad in the same way. The 
	# average transformation that align the left trajectory to the rigth trajectory (we use the notation in 
	#Horn 1987 and Picco 2015) is only true for small rotations and it is important to minimise inaccuracies
	#that can derive from the approximation of the rotation and traslation. For more details see 
	#Picco et al. 2015, Material and Methods, Two color alignment procedure, Estimate of the average trasformations).
	#################################################################################################################

	target_trajectory = Traj()
	target_trajectory.load( path_target )

	if ( fimax1 ) :

		print( 'fimax1 = True ; the software uses only the information of the target trajectory up to its peak of fluorescence intensity.' )
		
		t1 = target_trajectory.fimax( fimax_filter )
		t1.start( float( target_trajectory.annotations( 'raw_traj_starts_mean' ) ) )
	
	else :

		t1 = target_trajectory
		t1.start( float( target_trajectory.annotations( 'raw_traj_starts_mean' ) ) )
		t1.end( float( target_trajectory.annotations( 'raw_traj_ends_mean' ) ) )

	t1_center_mass = t1.center_mass()
	t1.translate( - t1_center_mass )

	#load the reference trajectory
	reference_trajectory = Traj()
	reference_trajectory.load( path_reference )
	
	if ( fimax2 ) :
		
		print( 'fimax2 = True ; the software uses only the information of the reference trajectory up to its peak of fluorescence intensity.' )
		
		t2 = reference_trajectory.fimax( fimax_filter )
		t2.start( float( target_trajectory.annotations( 'raw_traj_starts_mean' ) ) )
	
	else :

		t2 = reference_trajectory
		t2.start( float( target_trajectory.annotations( 'raw_traj_starts_mean' ) ) )
		t2.end( float( target_trajectory.annotations( 'raw_traj_ends_mean' ) ) )

	t2_center_mass = t2.center_mass()
	t2.translate( - t2_center_mass )
	
	l = len( ch1 )
	
	#control that the dataset of loaded trajectories is complete
	if l != len( ch2 ) : raise IndexError( 'The number of trajectories for ch1 and for ch2 differ.' )

	#define the dictionary where the transformations will be stored
	T = { 'angle' : [] , 'translation' : [] , 'lag' : [] }

	#compute the transformations that align t1 and t2 together.
	for i in range( l ) :

		print( "Align " + path_target + " to " + ch1[ i ].annotations( 'file' ) + " and " + path_reference + " to " + ch2[ i ].annotations( 'file' ) ) 

		#spline the target trajectories, to reduce the noise
		if ( fimax1 ) :
			spline_t1 , spline_ch1 = spline( t1 , ch1[ i ].fimax( fimax_filter ) )
		else :
			spline_t1 , spline_ch1 = spline( t1 , ch1[ i ] )

		#lag t1
		ch1_lag = cc( spline_t1 , spline_ch1 )
		spline_ch1.input_values( 't' , spline_ch1.t() + ch1_lag )
		
		#spline the reference trajectories, to reduce the noise
		if ( fimax2 ) :
			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ].fimax( fimax_filter ) )
		else :
			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ] )

		#lag t2
		ch2_lag = cc( spline_t2 , spline_ch2 )
		spline_ch2.input_values( 't' , spline_ch2.t() + ch2_lag )

		#unify the start and the end of the trajectory splines that are paired to compute the rotation and translation.
		unify_start_and_end( spline_t1 , spline_ch1 )
		unify_start_and_end( spline_t2 , spline_ch2 )
	
		#NOTE: the weight used in Picco et al., 2015 is slightly different. To use the same weight one should replace spline_t1.f() with spline_t1.f() / ( spline_t1.coord_err()[ 0 ] * spline_t1.coord_err()[ 1 ] )
		align_ch1_to_t1 = MSD( spline_t1 , spline_ch1 ) 
		align_ch2_to_t2 = MSD( spline_t2 , spline_ch2 )

		#The tranformation that aligns t1 to t2 will be the transformation that align ch2 to t2 and the 
		#inverse of the transformation that aligns ch1 to t1.
		#
		# R_2 @ R_1^{-1} @ ( t1 - t1.center_mass() ) + R_2 @ ( ch1.center_mass() - ch2.center_mass() ) + t2.center_mass()
		#
		#As the mean in MSD is weighted (see MSD in trajalign/average.py) the equation becomes
		#
		# R_2 @ R_1^{-1} @ ( t1 - align_ch1_to_t1[ 'rc' ] ) + R_2 @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] ) + align_ch2_to_t2[ 'rc' ] 
		#
		#where align_ch1_to_t1[ 'rc' ], align_ch1_to_t1[ 'lc' ], align_ch2_to_t2[ 'rc' ] and align_ch2_to_t2[ 'lc' ] are 
		#the estimates of the center of masses with the weight mean convention used in MSD.
		#Finally, the target and reference trajectory were initially shifted by 
		#
		# - t2_center_mass 
		#
		#Therefore, the final transformation that align the target trajectory to the reference trajectory must be corrected for this initial shifts
		#
		# R_2 @ R_1^{-1} @ ( t1 - align_ch1_to_t1[ 'rc' ] ) + R_2 @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] ) + align_ch2_to_t2[ 'rc' ] +
		# + t2_center_mass
		#
		#NOTE: in eLife we used the geometrical center of mass, t1.center_mass(), and not the 
		#approximation of the center of mass that best align t1 and ch1 under the weight convention in MSD, which is align_ch1_to_t1[ 'rc' ].
		#Therefore, in Picco et al, 2015
		#
		#	- R( T[ 'angle' ][ -1 ] ) @ align_ch1_to_t1[ 'rc' ]
		#
		#woud become 
		#
		#	- R( T[ 'angle' ][ -1 ] ) @ t1.center_mass()
		#
		
		#Compute the angle as the atan2 of the sin( align_ch2_to_t2[ 'angle' ] - align_ch1_to_t1[ 'angle' ] ) 
		#and cos( align_ch2_to_t2[ 'angle' ] - align_ch1_to_t1[ 'angle' ] ) 
		a = np.sin( align_ch2_to_t2[ 'angle' ] ) * np.cos( align_ch1_to_t1[ 'angle' ] ) -  np.cos( align_ch2_to_t2[ 'angle' ] ) * np.sin( align_ch1_to_t1[ 'angle' ] )  
		b = np.cos( align_ch2_to_t2[ 'angle' ] ) * np.cos( align_ch1_to_t1[ 'angle' ] ) +  np.sin( align_ch2_to_t2[ 'angle' ] ) * np.sin( align_ch1_to_t1[ 'angle' ] )  
		T[ 'angle' ].append( np.arctan2( a , b ) )
		T[ 'translation' ].append( np.array( 
				- R( T[ 'angle' ][ -1 ] ) @ align_ch1_to_t1[ 'rc' ]\
						+ R( align_ch2_to_t2[ 'angle' ] ) @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] )\
						+ align_ch2_to_t2[ 'rc' ] + t2_center_mass 
				)[ 0 ] ) #the [ 0 ] is because otherwise it would be [[ x , y ]] instead of [ x , y ]
		T[ 'lag' ].append( ch2_lag - ch1_lag )

	#compute the median and the standard error (SE) of the transformations.
	#NOTE that if fimax2 is used, the center of mass of reference trajectory does not 
	#correspond to the center of mass of the trajectory to which the target trajectory 
	#is aligned. The target trajectory, in fact, is aligned to the center of mass of
	T_median = { 
			'angle' : np.median( T[ 'angle' ] ) ,
			'angle_SE' : nanMAD( T[ 'angle' ] ) / np.sqrt( l ) ,
			'translation' : [
				np.median( [ T[ 'translation' ][ i ][ 0 ] for i in range( l ) ] ) ,
				np.median( [ T[ 'translation' ][ i ][ 1 ] for i in range( l ) ] )
				] ,
			'translation_SE' : [
				nanMAD( [ T[ 'translation' ][ i ][ 0 ] for i in range( l ) ] ) / np.sqrt( l ),
				nanMAD( [ T[ 'translation' ][ i ][ 1 ] for i in range( l ) ] ) / np.sqrt( l ) 
				] ,
			'lag' : np.median( T[ 'lag' ] ) , 
			'lag_SE' : nanMAD( T[ 'lag' ] ) / np.sqrt( l ) ,
			'n' : l
			}

	target_trajectory.rotate( T_median[ 'angle' ] , 
			angle_err = T_median[ 'angle_SE' ]
			)
	target_trajectory.translate( T_median[ 'translation' ] , 
			v_err = ( T_median[ 'translation_SE' ][ 0 ] , T_median[ 'translation_SE' ][ 1 ] )
			)
	target_trajectory.input_values( 't' , target_trajectory.t() + T_median[ 'lag' ] )

	dot_positions = [ i for i in range(len( path_target )) if path_target[i] == '.' ]
	file_ending = dot_positions[ len(dot_positions) - 1 ] #there could be more than one dot in the file name. Pick the last.
	file_name =  path_target[ 0 : file_ending ] + '_aligned' + path_target[ file_ending : len( path_target ) ]

	# annotations
	target_trajectory.annotations( 'aligned_to' , str( path_reference ) )
	target_trajectory.annotations( 'original_file' , str( path_target ) )
	target_trajectory.annotations( 'alignment_angle' , str( T_median[ 'angle' ] ) + ' rad' )
	target_trajectory.annotations( 'alignment_angle_SE' , str( T_median[ 'angle_SE' ] ) + ' rad' )
	target_trajectory.annotations( 'alignment_translation' , str( T_median[ 'translation' ] ) + ' ' + target_trajectory.annotations()[ 'coord_unit' ] )
	target_trajectory.annotations( 'alignment_translation_SE' , str( T_median[ 'translation_SE' ] ) + ' ' + target_trajectory.annotations()[ 'coord_unit' ] )
	target_trajectory.annotations( 'alignment_lag' , str( T_median[ 'lag' ] ) + ' ' + target_trajectory.annotations()[ 't_unit' ] )
	target_trajectory.annotations( 'alignment_lag_SE' , str( T_median[ 'lag_SE' ] ) + ' ' + target_trajectory.annotations()[ 't_unit' ] )

	target_trajectory.save( file_name )

	print( 'The trajectory aligned to ' + path_reference + ' has been saved as ' + file_name )

	#-------------------------END-OF-align-DEFINITION--------------------------------

def average_ch1( path_reference , ch1 , ch2 , output_file = 'average' , median = False , unify_start_end = False , max_frame = [] ,  fimax = False , fimax_filter = [ -3/35 , 12/35 , 17/35 , 12/35 , -3/35 ] ):

	"""
	average_one_channel( path_reference , ch1 , ch2 , fimax = False , fimax_filter ) : averages the trajectories listed in ch1. These trajectories  were acquired symultaneously to the trajectories in ch2, which are aligned to the reference trajectory identified by path_reference. The transformation that aligns the trajectories in ch2 to path reference is used to align the ch1 trajectories together and to compute then their average. fimax allows the user to use only the trajectory information up to the peak of flurescence intensity to compute the transformation.
	"""

	if unify_start_end and not max_frame :
	
		raise TypeError('You need to specify the max_frame if you want to unify the start and end')

	header()

	aligned_ch1 = []

	reference_trajectory = Traj()
	reference_trajectory.load( path_reference )

	if ( fimax ) :
		
		print( 'fimax = True ; the software uses only the information of the reference trajectory up to its peak of fluorescence intensity.' )
		
		t2 = reference_trajectory.fimax( fimax_filter )
	
	else :

		t2 = reference_trajectory

	t2_center_mass = t2.center_mass()
	t2.translate( - t2_center_mass )
	
	l = len( ch1 )
		
	#control that the dataset of loaded trajectories is complete
	if l != len( ch2 ) : raise IndexError( 'The number of trajectories for ch1 and for ch2 differ.' )

	#compute the transformations that align t1 and t2 together.
	trajectories_time_span = { 'old_start' : [], 'new_start' : [], 'old_end' : [], 'new_end' : []}
	for i in range( l ) :
		
		print( "Align " + ch1[ i ].annotations( 'file' ) + " and " + ch2[ i ].annotations( 'file' ) + " to " + path_reference ) 

		if ( fimax ) :

			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ].fimax( fimax_filter ) )

		else :

			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ] )

		#lag t2
		ch2_lag = cc( spline_t2 , spline_ch2 )
		spline_ch2.input_values( 't' , spline_ch2.t() + ch2_lag )

		#unify the start and the end of the trajectory splines that are paired to compute the rotation and translation.
		unify_start_and_end( spline_t2 , spline_ch2 )
		
		#compute the elements of the transformation that aligns ch2 to t2
		align_ch2_to_t2 = MSD( spline_t2 , spline_ch2 )

		aligned_ch1.append( cp.deepcopy( ch1[ i ] ) )

		aligned_ch1[ i ].rotate( align_ch2_to_t2[ 'angle' ])
		aligned_ch1[ i ].translate( np.array( align_ch2_to_t2[ 'rc' ] 
			- R( align_ch2_to_t2[ 'angle' ] ) @ ( align_ch2_to_t2[ 'lc' ] ) 
			+ t2_center_mass )[ 0 ] ) #the [ 0 ] is because otherwise it would be [[ x , y ]] instead of [ x , y ]
		aligned_ch1[ i ].input_values( 't' , aligned_ch1[ i ].t() + ch2_lag )

		trajectories_time_span[ 'old_start' ].append( ch1[ i ].start() )
		trajectories_time_span[ 'old_end' ].append( ch1[ i ].end() )
		trajectories_time_span[ 'new_start' ].append( aligned_ch1[ i ].start() )
		trajectories_time_span[ 'new_end' ].append( aligned_ch1[ i ].end() )
	
	#average the ch1 trajectories that have been aligned together

	if unify_start_end :

		mean_start , mean_end = compute_average_start_and_end( trajectories_time_span , aligned_ch1 , max_frame )
		
		#uniform start and end of aligned trajectories to mean_start and mean_end
		for j in range( l ):
		
			aligned_ch1[ j ].start( mean_start )
			aligned_ch1[ j ].end( mean_end )
	else :

		for j in range( l ):
		
			aligned_ch1[ j ].start( min( trajectories_time_span[ 'new_start' ] ) )
			aligned_ch1[ j ].end( max( trajectories_time_span[ 'new_end' ] ) )
				
	average_trajectory =  trajectory_average( aligned_ch1 , 0 , median , fimax )
	
	if not unify_start_end :

		average_trajectory.annotations( 'raw_traj_starts' , trajectories_time_span[ 'new_start' ] )
		average_trajectory.annotations( 'raw_traj_ends' , trajectories_time_span[ 'new_end' ] )
		average_trajectory.annotations( 'raw_traj_starts_mean' , np.mean( trajectories_time_span[ 'new_start' ] ) ) 
		average_trajectory.annotations( 'raw_traj_starts_std' , np.std( trajectories_time_span[ 'new_start' ] ) ) 
		average_trajectory.annotations( 'raw_traj_ends_mean' , np.mean( trajectories_time_span[ 'new_end' ] ) ) 
		average_trajectory.annotations( 'raw_traj_ends_std' , np.std( trajectories_time_span[ 'new_end' ] ) ) 

	#save the aligned_ch1 trajectories
	for i in range( l ) :
		filename = "./" + output_file + "/" + aligned_ch1[ i ].annotations( 'file' )
		if i == 0 :
			directory = os.path.dirname( filename )
			if not os.path.exists( directory ) :
				os.makedirs( directory )
		aligned_ch1[ i ].save( filename )
	
	mean_precision =  np.sqrt(
			np.nanmean( 
				average_trajectory.coord_err()[ 0 ] ** 2 + average_trajectory.coord_err()[ 1 ] ** 2 
				)
			)

	print('ALIGNMENT PRECISION:')
	print( mean_precision )

	average_trajectory.save( output_file )

