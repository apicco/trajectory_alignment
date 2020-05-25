import os 
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from trajalign.traj import Traj
from scipy.stats import  t as ttest
from scipy.stats import  f

import rpy2.robjects as r
from rpy2.robjects.packages import importr

MASS = importr( 'MASS' )
base = importr( 'base' )

def split_pt( path_input , path_outputs , i0 = -1 , pattern = '%% Trajectory' ) :

	i = i0

	with open( path_input , 'r' ) as f :

		for line in f :
			
			if pattern in line :

				if ( i >= 0 ) & ( i0 == -1 ) : #then there is already a g open (see following lines). i0 not == -1 if there are already some files and the indexing has to start higher.

					g.close()

				i = i + 1

				i0 = -1 #i0 is restored to be -1, if it is not, so that from next iteration the g will be closed

			if i >= 0 :

				with open( path_outputs + '/trajectory_%06d' % i + '.txt' , 'a' ) as g :

					g.write( line ) 

	return i 
					

#def load_traj( path , pattern = 'x' , comment_char = '%' , **kwargs ) :
#	"""
#	load all the trajectories identified by pattern in path
#	"""
#
#	output = []
#
#	#list the trajectory files identified by pattern
#	r = [ f for f in os.listdir( path ) if pattern in f ]
#
#	#load the trajectory as a Traj object
#	for i in range( len( r ) ) :
#	
#		t = Traj()
#		t.load( path + '/' + r[ i ] , comment_char = '%' , **kwargs )
#		#t.extract( t.f() == t.f() ) #remove NA
#		output.append( t )
#
#	return output

def distance( a , b ) :

	return( np.sqrt( ( a[ 0 ] - b[ 0 ] ) ** 2 + ( a[ 1 ] - b[ 1 ] ) ** 2 ) )

def ecc( t ) :

	# Computes the eccentricity of the trajectory t

	u02 = t.u02()
	u20 = t.u20()
	u11 = t.u11()

	N = len( u02 )

	a = [ ( u02[i] + u20[i] ) / 2 for i in range( N ) ]
	b = [ np.sqrt( (u20[i] - u02[i]) ** 2 + 4 * u11[i] ** 2 ) / 2 for i in range( N ) ]
	
	# eccentricity, ellipse radii:
	l_1 = [ np.sqrt( a[i] + b[i] ) for i in range( N ) ]
	l_2 = [ np.sqrt( a[i] - b[i] ) for i in range( N ) ]
	
	# ration between ellipse radii:
	l_r = [ l_2[ i ] / l_1[ i ] for i in range( N ) ]
	# eccentricity: (because of their above definitions l_2[ i ] < l_1[ i ] for each i)
	e = [ np.sqrt( 1 - l_r[ i ] ) for i in range( N ) ] 
	
	return e  , [ l_1 , l_2 ]

def mean_centroid( x ) :

	return( [ np.nanmean( x.coord()[ 0 ] ) , np.nanmean( x.coord()[ 1 ] ) ] )

def eccStats( t , rt , m0 = 1 , c0 = 0 , maxit = 20 ) :

	# compare if the eccentricity values compute in the region where the spot
	# is quantified and in the surrounding region are falling on a line close to the 
	# diagonal: y = m0 * x + c0. The H0 is that y = m0 * x + c0 is sufficient to describe 
	# the correlation between the eccentricities (i.e. the spot needs to be kept).
	x , _ = ecc( t )
	y , _ = ecc( rt )

	# remove possible nan and select only eccentricities smaller than the median
	xx = [ x[ i ] for i in range( len( x ) ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) & ( x[ i ] <= np.nanmedian( x ) ) ) ]
	yy = [ y[ i ] for i in range( len( y ) ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) & ( x[ i ] <= np.nanmedian( x ) ) ) ]

	chisq = sum( [ ( yy[ i ] - xx[ i ] ) **2 / xx[ i ] for i in range( len( xx ) ) ] )

	return chisq

def ichose( tt , rtt , image_shape, image_len, pval_m = 0.1 , pval_c = 0.1 , pval_F = 1 , maxit = 100 , d0 = 10 , t0 = 0 ) :
	"""
	ichose( tt , rtt , image_shape, image_len, pval_m = 0.1 , pval_c = 0.1 , pval_F = 1 , maxit = 100 , d0 = 10 , t0 = 0 )
	select the trajectories in the trajectory list 'tt' whose spots have eccentricities that do not change if measured in
	a larger region around the spot. These eccentricity values are inputed through the trajectory list 'rtt'. 
	- image_shape, define the shape of the image so that d0 (see below) can be used
	- image_len, define the length of the image so that t0 (see below) can be used
	- pval_m and pval_c define the cutoffs on the t-tests on the interpolation of the eccentricities. The higher the pvalues, the 
	more stringent the spot selection
	- pval_F, deprecated
	- maxit, max number of iterations in rlm (R)
	- d0, the region on the image border where trajectories are rejected by default (too close to the image edges)
	- t0, trajectories starting before frame t0 and ending after frame image_len - t0 - 1 are rejected
	"""

	# d0 sets the minimal distance from the image border
	output_tt = []
	output_rtt = []

	l = len( tt ) 

	if not l == len( rtt ) :

		raise TypeError( 'the list of trajectories tt and rtt have different lengths' )

	for i in range( l ) :

		pm , m , pc , c , pf = eccStats( tt[ i ] , rtt[ i ] , maxit = maxit )

		"The H0 is that the data are well described by "
		if ( ( pm > pval_m ) & ( pc > pval_c ) & ( pf < pval_F ) ) : 

			mc = mean_centroid( tt[ i ] )

			if ( ( mc[ 0 ] > d0 ) & ( mc[ 1 ] > d0 ) ) :
				
				if ( ( mc[ 0 ] < ( image_shape[ 0 ] -  d0 ) ) & ( mc[ 1 ] < ( image_shape[ 1 ] - d0 ) ) ) :

					if ( ( tt[ i ].frames( 0 ) > t0 ) & ( tt[ i ].frames( len( tt[ i ] ) - 1 ) < image_len - 1 - t0 ) ) :
					
						# annotations useful to check selection parameters
						tt[ i ].annotations( 'eccentricity_m' , str( m[ 0 ] ) )
						tt[ i ].annotations( 'eccentricity_c' , str( c[ 0 ] ) )
						
						tt[ i ].annotations( 'eccentricity_pval_m' , str( pm ) )
						tt[ i ].annotations( 'eccentricity_pval_c' , str( pc ) )
						tt[ i ].annotations( 'eccentricity_pval_F' , str( pf ) )
	
						tt[ i ].annotations( 'threshold_pval_m' , str( pval_m ) )
						tt[ i ].annotations( 'threshold_pval_c' , str( pval_c ) )
						tt[ i ].annotations( 'threshold_pval_F' , str( pval_F ) )
						
						output_tt.append( tt[ i ] )
						output_rtt.append( rtt[ i ] )

	return output_tt , output_rtt

def save_directory( tt , directory_path ) :
	"""
	save_directory( tt , directory_path ) saves the trajectories in the list 'tt' as files in the directory path named 'directory_path'. The file name and extension will be the same as in the .annotations()[ 'file' ]
	"""
	

	if not os.path.exists( directory_path ) :
		
		os.makedirs( directory_path )

	for t in tt :

		file_name = t.annotations()[ 'file' ]
		t.save( directory_path + '/' + file_name )


def plot_traj( tt , f , what , ms = 30 , lw = 3 ) :

	cmap = cm.get_cmap( 'prism' , len( tt ) ) #color map
	
	shift = 0#0.5 #correct PT shift

	l = len( tt )

	xlim_max = []
	ylim_min = []
	ylim_max = []

	for j in range( l ) :
	
		c = cmap( j / ( l ) ) #colormap
		t = tt[ j ]
		

		if ( f + 1 >= t.frames()[ 0 ] ) & ( f + 1 <= t.frames()[ -1 ] ) :

			#selected frames
#			s = list( range( t.frames()[ 0 ] , f + 1 ) ) 
#			sel = [ i for i in range( len( s ) ) if s[ i ] in t.frames() ]
	
			sel = [ i for i in range( len( t ) ) if t.frames()[ i ] <= f + 1 ]
	
			u = t.extract( sel )

			if ( what == 'coord' ) :
		
				plt.plot( u.coord()[ 1 ][ -1 ] + shift , u.coord()[ 0 ][ -1 ] - shift , 'o' , color = c , markersize = ms , markerfacecolor = 'none' )
				plt.plot( u.coord()[ 1 ] + shift, u.coord()[ 0 ] - shift , '-' , color = c , linewidth = lw )

			else :
		
				print( 'something' ) 
				#es = eccStats( t , rt )
				#plt.plot( es[ 0 ] , es[ 1 ] , 'o' , color = c , markersize = ms )

