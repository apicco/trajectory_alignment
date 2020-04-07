import os 
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from trajalign.traj import Traj
from skimage.external import tifffile as tiff
from scipy.stats import f , norm


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

	rs1 = [ ( y[ i ] - np.nanmean( y ) ) ** 2 for i in range( n ) if y[ i ] == y[ i ] ]
	rss1 = sum( rs1 )
	rs2 = [ ( y[ i ] - x[ i ] ) ** 2 for i in range( n ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] # predicted y is x, because the model function is y = x 
	rss2 = sum( rs2 )

	N = len( rs1 )

	if N == len( rs2 ) : 

		if ( ( p2 - p1 ) > 0 ) & ( ( N - p2 ) > 0 ) :
		
			F = ( ( rss1 - rss2 ) / ( p2 - p1 ) ) / ( rss2 / ( N - p2 ) )
	
			p = 1 - f.cdf( F , p2 - p1 , N - p2 )
	
		else :
	
			p = F = np.nan
	
		return p , F #H0: model does not provide a better description of the data than the restricted model (i.e. average, where all possible parameters are restricted to 0)

	else :

		raise TypeError( 'There is an error in your data, the number of NaN differs between y and x' )

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

	# d0 sets the minimal distance from the image border

	output_tt = []
	output_rtt = []

	l = len( tt ) 

	if not l == len( rtt ) :

		raise TypeError( 'the list of trajectories tt and rtt have different lengths' )

	for i in range( l ) :

		p = eccStats( tt[ i ] , rtt[ i ] )
	
		if p <= pval : 

			m = mean_centroid( tt[ i ] )

			if ( ( m[ 0 ] > d0 ) & ( m[ 1 ] > d0 ) ) :
				
				if ( ( m[ 0 ] < ( image_shape[ 0 ] -  d0 ) ) & ( m[ 1 ] < ( image_shape[ 1 ] - d0 ) ) ) :

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
	print( 'mean : ' + str( np.nanmean( e ) ) )
	print( 'std : ' + str( np.nanstd( e ) ) )
	print( 'p : ' + str( p ) )
	
	return l_1 , l_2 , e , p 


