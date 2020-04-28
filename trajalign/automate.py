import os 
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from trajalign.traj import Traj
from skimage.external import tifffile as tiff
from scipy.stats import  t as ttest
from scipy.stats import norm , t , f 


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

def eccStats( t , rt , m0 = 1 , c0 = 0 ) :

	# compare if the eccentricity values compute in the region where the spot
	# is quantified and in the surrounding region are falling on a line close to the 
	# diagonal: y = m0 * x + c0. The H0 is that y = m0 * x + c0 is sufficient to describe 
	# the correlation between the eccentricities (i.e. the spot needs to be kept).
	x , _ = ecc( t )
	y , _ = ecc( rt )
	
	# remove possible nan
	xx = [ x[ i ] for i in range( len( x ) ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) & ( x[ i ] <= np.nanmedian( x ) ) ) ]
	yy = [ y[ i ] for i in range( len( y ) ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) & ( y[ i ] <= np.nanmedian( y ) ) ) ]

	# degree of freedom are n - 2 (m + c, two parameters to be fixed)
	n = len( xx )
	df = n - 2 

	# linear regression of order 1, output covariance matrix whose diagonal elements are
	# the variance used to compute the SE over the estimates of m and c, which are then
	# used to compute a t test
	if n > 4 : #he number of data points must exceed order + 2 ( + 1? )
		
		fit = np.polyfit( xx , yy , 1 , cov = True )
	
		# SE are computed according to the definition used in python with the divident (n - 2)
		# in the hope that small trajectories are penalized
		c = fit[ 0 ][ 1 ]
		sc = np.sqrt( fit[ 1 ][ 1 , 1 ] )
		m = fit[ 0 ][ 0 ]
		sm = np.sqrt( fit[ 1 ][ 0 , 0 ] )
	
		# the t variables are
		t_c = np.abs( ( c - c0 ) / sc )
		t_m = np.abs( ( m - m0 ) / sm )
	
		# the p values, testing the H0 that the data are sufficiently described by the simple model defined
		# by the parameters m0 and c0 rather than by the one from the interpolation.
		p_c = 1 - ttest.cdf( t_c , df ) + ttest.cdf( -t_c , df )
		p_m = 1 - ttest.cdf( t_m , df ) + ttest.cdf( -t_m , df )

	else : 

		p_m = p_c = m = sm = c = sc = np.nan
	
	# Also, perform an F-test to exclude that data are not correlated, with:
	p1 = 1
	p2 = 2

	rs1 = [ ( yy[ i ] - np.mean( xx[ i ] ) ) ** 2 for i in range( n ) ]
	rss1 = sum( rs1 )
	rs2 = [ ( yy[ i ] - m * xx[ i ] - c ) ** 2 for i in range( n ) ] # predicted y is x, because the model function is y = x 
	rss2 = sum( rs2 )

	if   n > 2 :
	
		F = ( ( rss1 - rss2 ) / ( p2 - p1 ) ) / ( rss2 / ( n - p2 ) )
		pf = 1 - f.cdf( F , p2 - p1 , n - p2 )

	else :

		pf = F = np.nan

	return p_m , [ m , sm ] , p_c , [ c , sc ] , pf

def ichose( tt , rtt , image_shape, pval_m = 0.1 , pval_c = 0.0 , pval_F = 0.0001 , d0 = 10 ) :

	# d0 sets the minimal distance from the image border
	output_tt = []
	output_rtt = []

	l = len( tt ) 

	if not l == len( rtt ) :

		raise TypeError( 'the list of trajectories tt and rtt have different lengths' )

	for i in range( l ) :

		pm , m , pc , c , pf = eccStats( tt[ i ] , rtt[ i ] )

		"The H0 is that the data are well described by "
		if ( ( pm > pval_m ) & ( pc > pval_c ) & ( pf < pval_F ) ) : 

			mc = mean_centroid( tt[ i ] )

			if ( ( mc[ 0 ] > d0 ) & ( mc[ 1 ] > d0 ) ) :
				
				if ( ( mc[ 0 ] < ( image_shape[ 0 ] -  d0 ) ) & ( mc[ 1 ] < ( image_shape[ 1 ] - d0 ) ) ) :
				
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

def tmp_ecc( t ) :

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
	
	N = len( l_1 )
	# ration between ellipse radii:
	l_r = [ l_2[ i ] / l_1[ i ] for i in range( N ) ]
	# eccentricity: (because of their above definitions l_2[ i ] < l_1[ i ] for each i)
	e = [ np.sqrt( 1 - l_r[ i ] ) for i in range( N ) ] 

	p =  norm.cdf( 0 , np.nanmean( e ) , np.nanstd( e ) )
	
	return l_1 , l_2 , e , p 


