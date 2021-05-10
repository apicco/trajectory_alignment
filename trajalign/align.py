# All the software here is distributed under the terms of the GNU General Public License Version 3, June 2007. 
# Trajalign is a free software and comes with ABSOLUTELY NO WARRANTY.
# 
# You are welcome to redistribute the software. However, we appreciate is use of such software would result in citations of 
# Picco, A., Kaksonen, M., _Precise tracking of the dynamics of multiple proteins in endocytic events_,  Methods in Cell Biology, Vol. 139, pages 51-68 (2017)
# http://www.sciencedirect.com/science/article/pii/S0091679X16301546
# 
# Author: Andrea Picco (https://github.com/apicco)
# Year: 2017

from trajalign.traj import Traj
from trajalign.average import load_directory
from trajalign.average import MSD
from trajalign.average import nanMAD 
from trajalign.average import header
from trajalign.average import unified_start , unified_end
from scipy.interpolate import UnivariateSpline #need to install py35-scikit-learn
import numpy as np
import copy as cp

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
				unit = to_interpolate.annotations()[ 't_unit' ]
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

def align( path_target , path_reference , ch1 , ch2 , fimax1 = False , fimax2 = False , fimax_filter = [ -3/35 , 12/35 , 17/35 , 12/35 , -3/35 ] , unify_start_end_in_alignment = True , unify_start_end_in_output = False ):

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

	header() 

	target_trajectory = Traj()
	target_trajectory.load( path_target )

	reference_trajectory = Traj()
	reference_trajectory.load( path_reference )
	
	#################################################################################################################
	#average trajectories are centered on their center of mass and must have been previously lied down 
	#(lie_down function in trajalign/average.py) so that they are orientad in the same way. The 
	# average transformation that align the left trajectory to the rigth trajectory (we use the notation in 
	#Horn 1987 and Picco 2015) is only true for small rotations and it is important to minimise inaccuracies
	#that can derive from the approximation of the rotation and traslation. For more details see 
	#Picco et al. 2015, Material and Methods, Two color alignment procedure, Estimate of the average trasformations).
	#################################################################################################################
	
	if ( fimax1 ) :

		print( 'fimax1 = True ; the software uses only the information of the target trajectory up to its peak of fluorescence intensity.' )
		
		t1 = target_trajectory.fimax( fimax_filter )
	
	else :

		t1 = cp.deepcopy( target_trajectory )

	if ( fimax2 ) :
		
		print( 'fimax2 = True ; the software uses only the information of the reference trajectory up to its peak of fluorescence intensity.' )
		
		t2 = reference_trajectory.fimax( fimax_filter )
	
	else :

		t2 = cp.deepcopy( reference_trajectory )

	print( "unify_start_end_in_output : " + str( unify_start_end_in_output ) )
	
	if unify_start_end_in_alignment :

		#if the input average trajectories have not unified start and end (unify_start_end = False), then unify the start and end for the sake of the alignment.

		if t1.annotations()[ 'unify_start_end' ] is 'False' :
			
			t1.start( unified_start( t1 ) )
			t1.end( unified_end( t1 ) )

			print( '\nunify_start_end_in_alignment = True ; the target average trajectory new start and end are ' + str( unify_start( t1 ) ) + ' ' + t1.annotations()[ 't_unit' ] + ' and ' + str( unify_end( t1 ) ) + ' ' + t1.annotations()[ 't_unit' ] + '.\n' )

		else :

			print( '\nunify_start_end_in_alignment = ' + str( unify_start_end_in_alignment ) + '. The target average trajectory had already unified start and end values.\n' )
	
		
		if t2.annotations()[ 'unify_start_end' ] is 'False' :
			
			t2.start( unified_start( t2 ) )
			t2.end( unified_end( t2 ) )

			print( 'unify_start_end_in_alignment = True ; the reference average trajectory new start and end are ' + str( unify_start( t2 ) ) + ' ' + t2.annotations()[ 't_unit' ] + ' and ' + str( unify_end( t2 ) ) + ' ' + t2.annotations()[ 't_unit' ] + '.\n' )
		
		else :

			print( 'unify_start_end_in_alignment = ' + str( unify_start_end_in_alignment ) + '. The reference average trajectory had already unified start and end values.\n' )
	
	else :	
	
		if t1.annotations()[ 'unify_start_end' ] is 'True' :

			print( '\nunify_start_end_in_alignment = False but the target average trajectory was computed with unify_start_end = True. I cannot perform this operation. unify_start_end_in_alignment set to True for the target trajectory.\n' )
		
		else :

			print( '\nunify_start_end_in_alignment = ' + str( unify_start_end_in_alignment ) + '\n' )
	
		if t2.annotations()[ 'unify_start_end' ] is 'True' :

			print( 'unify_start_end_in_alignment = False but the reference average trajectory was computed with unify_start_end = True. I cannot perform this operation. unify_start_end_in_alignment set to True for the reference trajectory\n' )
		
		else :

			print( 'unify_start_end_in_alignment = ' + str( unify_start_end_in_alignment ) + '\n' )

	t1_center_mass = t1.center_mass()
	t1.translate( - t1_center_mass )

	t2_center_mass = t2.center_mass()
	t2.translate( - t2_center_mass )
	
	print( "------------------DEBUG---------------------------")
	print( "t1 cm = " + str( t1_center_mass ) + "; target_trajectory cm =" + str( target_trajectory.center_mass() ) )
	print( "t2 cm = " + str( t2_center_mass ) + "; reference_trajectory cm =" + str( reference_trajectory.center_mass() ) )
	print( "------------------DEBUG---------------------------")

	l = len( ch1 )
	
	#control that the dataset of loaded trajectories is complete
	if l != len( ch2 ) : raise IndexError( 'The number of trajectories for ch1 and for ch2 differ.' )

	#define the dictionary where the transformations will be stored
	T = { 'angle' : [] , 'translation' : [] , 'lag' : [] }

	#compute the transformations that align t1 and t2 together.
	for i in range( l ) :

		print( "Align " + path_target + " to " + ch1[ i ].annotations()[ 'file' ] + " and " + path_reference + " to " + ch2[ i ].annotations()[ 'file' ] ) 

		#spline the trajectories, to reduce the noise
		if ( fimax1 ) :
			spline_t1 , spline_ch1 = spline( t1 , ch1[ i ].fimax( fimax_filter ) )
		else :
			spline_t1 , spline_ch1 = spline( t1 , ch1[ i ] )

		if ( fimax2 ) :
			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ].fimax( fimax_filter ) )
		else :
			spline_t2 , spline_ch2 = spline( t2 , ch2[ i ] )

		#lag t1
		ch1_lag = cc( spline_t1 , spline_ch1 )
		spline_ch1.input_values( 't' , spline_ch1.t() + ch1_lag )
		
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
		#
		# - t2_center_mass 
		#
		#and
		#
		# - t1_center_mass 
		#
		#Therefore, the final transformation that align the target trajectory to the reference trajectory must be corrected for this initial shifts
		#
		# R_2 @ R_1^{-1} @ ( t1 - align_ch1_to_t1[ 'rc' ] ) + R_2 @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] ) + align_ch2_to_t2[ 'rc' ] +
		# + t2_center_mass + t1_center_mass
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
						+ align_ch2_to_t2[ 'rc' ] #+ t2_center_mass + t1_center_mass
				)[ 0 ] ) #the [ 0 ] is because otherwise it would be [[ x , y ]] instead of [ x , y ]
#bkp		T[ 'translation' ].append( np.array( 
#bkp				- R( T[ 'angle' ][ -1 ] ) @ align_ch1_to_t1[ 'rc' ]\
#bkp						+ R( align_ch2_to_t2[ 'angle' ] ) @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] )\
#bkp						+ align_ch2_to_t2[ 'rc' ] #+ t2_center_mass + t1_center_mass
#bkp				)[ 0 ] ) #the [ 0 ] is because otherwise it would be [[ x , y ]] instead of [ x , y ]
		T[ 'lag' ].append( ch2_lag - ch1_lag )

		#debug
		print( "------------------DEBUG---------------------------")
		tmp1 =  np.array( 
				- R( T[ 'angle' ][ -1 ] ) @ align_ch1_to_t1[ 'rc' ]\
						+ R( align_ch2_to_t2[ 'angle' ] ) @ ( align_ch1_to_t1[ 'lc' ] - align_ch2_to_t2[ 'lc' ] )\
						+ align_ch2_to_t2[ 'rc' ] 
				)[ 0 ] #the [ 0 ] is because otherwise it would be [[ x , y ]] instead of [ x , y ]
		print( "T" )
		print( T[ 'translation' ][ len(  T[ 'translation' ] ) - 1 ] )
		print( "T tmp1" )
		print( tmp1 )

		print( "center mass t1: "+ str( t1_center_mass ) ) 
		print( "target center mass : "+ str( target_trajectory.center_mass() ) )
		print( "center mass t2: "+ str( t2_center_mass ) ) 
		print( "ref center mass : "+ str( reference_trajectory.center_mass() ) )
		print( "------------------DEBUG---------------------------")
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
#	target_trajectory.annotations( 'starting_center_mass' , str( t1_center_mass )  + ' ' + target_trajectory.annotations()[ 'coord_unit' ] )
	target_trajectory.annotations( 'alignment_lag' , str( T_median[ 'lag' ] ) + ' ' + target_trajectory.annotations()[ 't_unit' ] )
	target_trajectory.annotations( 'alignment_lag_SE' , str( T_median[ 'lag_SE' ] ) + ' ' + target_trajectory.annotations()[ 't_unit' ] )
	target_trajectory.annotations( 'unify_start_end_in_alignment_output' , str( unify_start_end_in_output ) )

	# update the mean_starts and mean_ends, which are used for unify_start and unify_end.
	target_trajectory.annotations( 'mean_starts' , str( float( target_trajectory.annotations()[ 'mean_starts' ] ) + T_median[ 'lag' ] ) )
	target_trajectory.annotations( 'mean_ends' , str( float( target_trajectory.annotations()[ 'mean_ends' ] ) + T_median[ 'lag' ] ) )

	if unify_start_end_in_output :

		target_trajectory.start( unified_start( target_trajectory ) )
		target_trajectory.end( unified_end( target_trajectory ) )

	target_trajectory.save( file_name )

	print( 'The trajectory aligned to ' + path_reference + ' has been saved as ' + file_name )

