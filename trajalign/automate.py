import os 
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from trajalign.traj import Traj
from skimage.external import tifffile as tiff
from scipy.stats import f


def split_pt( path_input , path_outputs , i0 = -1 , pattern = '%% Trajectory' ) :

	i = i0

	with open( path_input , 'r' ) as f :

		for line in f :
			
			if pattern in line :

				if ( i >= 0 ) & ( i0 == -1 ) : #then there is already a g open. i0 not == -1 if there are already some files and the indexing has to start higher.

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

	# eccentricity, ellipse radii:
	l_1 = [ np.sqrt( ( u02[i] + u20[i] + np.sqrt( (u20[i] - u02[i]) ** 2 + 4 * u11[i] ** 2 ) ) / 2 ) for i in range( len( u02 ) ) ]
	l_2 = [ np.sqrt( ( u02[i] + u20[i] - np.sqrt( (u20[i] - u02[i]) ** 2 + 4 * u11[i] ** 2 ) ) / 2 ) for i in range( len( u02 ) ) ]
	
	N = len( l_1 )
	# ration between ellipse radii:
	l_r = [ l_2[ i ] / l_1[ i ] for i in range( N ) ]
	# eccentricity: (because of their above definitions l_2[ i ] < l_1[ i ] for each i)
	e = [ np.sqrt( 1 - l_r[ i ] ) for i in range( N ) ] 
	
	return e 

def yxF( x , y , p1 = 1 , p2 = 2 ) : 

	n = len( x )
	
	if not n == len( y ) :

		raise AttributeError( 'lenght of x and y must be equal' ) 

	rss1 = sum( [ ( y[ i ] - np.nanmean( y ) ) ** 2 for i in range( n ) if y[ i ] == y[ i ] ] )
	rss2 = sum( [ ( y[ i ] - x[ i ] ) ** 2 for i in range( n ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] ) # predicted y is x, because the model function is y = x 

	if ( ( p2 - p1 ) > 0 ) & ( ( n - p2 ) > 0 ) :
	
		F = ( ( rss1 - rss2 ) / ( p2 - p1 ) ) / ( rss2 / ( n - p2 ) )

		p = 1 - f.cdf( F , p2 - p1 , n - p2 )

	else :

		p = F = np.nan

	return p , F

def mean_centroid( x ) :

	return( [ np.nanmean( x.coord()[ 0 ] ) , np.nanmean( x.coord()[ 1 ] ) ] )

def eccStats( t , rt ) :
	
	x = ecc( t )
	y = ecc( rt )

	p , F = yxF( x , y )

	return p
	#TO-TRY: try to use shannon entropy with residuals or ratios. both residuals
	#and ratios should approx same values (i.e. shannon entropy max) in distribution
	# residuals from y = x
	##R = [ ( y[ i ] - x[ i ] ) ** 2 for i in range( len( x ) ) ] 
	##return  [ R , [ np.log( p ) ] * len( R ) ] #[ np.median( r ) , np.median( R ) ]

def ichose( tt , rtt , image_shape, pval = 0.01 , d0 = 10 ) :

	output = []

	l = len( tt ) 

	if not l == len( rtt ) :

		raise TypeError( 'the list of trajectories tt and rtt have different lengths' )

	for i in range( l ) :

		p = eccStats( tt[ i ] , rtt[ i ] )
	
		if p <= pval : 

			m = mean_centroid( tt[ i ] )

			if ( ( m[ 0 ] > d0 ) & ( m[ 1 ] > d0 ) ) :
				
				if ( ( m[ 0 ] < ( image_shape[ 0 ] -  d0 ) ) & ( m[ 1 ] < ( image_shape[ 1 ] - d0 ) ) ) :

					output.append( tt[ i ] )

	return output

def save_directory( tt , directory_path ) :
	"""
	save_directory( tt , directory_path ) saves the trajectories in the list 'tt' as files in the directory path named 'directory_path'. The file name and extension will be the same as in the .annotations()[ 'file' ]
	"""
	
	d = os.path.dirname( directory_path )

	if not os.path.exists( d ) :
		
		os.makedirs( d )

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

#	if not ( what == 'coord' ) : 
	
#		plt.plot( ( 0 , 0.5 ) , ( np.log( 0.05 ) , np.log( 0.05 ) ) , ls = '--' , lw = 1 , color = 'red' )
#		plt.plot( ( 0 , 0.5 ) , ( np.log( 0.01 ) , np.log( 0.01 ) ) , ls = '--' , lw = 1 , color = 'green' )
#		plt.plot( ( 0 , 0.5 ) , ( np.log( 0.001 ) , np.log( 0.001 ) ) , ls = '--' , lw = 1 , color = 'blue' )
#		plt.plot( ( 0 , 0.5 ) , ( np.log( 0.0001 ) , np.log( 0.0001 ) ) , ls = '--' , lw = 1 , color = 'black' )
#
#		plt.xlabel( 'Square residuals' )
#		plt.ylabel( '$ln( p )$' )
