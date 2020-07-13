import os 
import warnings 
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from trajalign.traj import Traj
from scipy.stats import chi2 , linregress

import rpy2.robjects as r
from rpy2.robjects.packages import importr

MASS = importr( 'MASS' )
base = importr( 'base' )

def split_pt( path_input , path_outputs , i0 = 0 , pattern = '%% Trajectory' ) :

	i = i0
	writing = False 

	with open( path_input , 'r' ) as f :

		for line in f :

			if pattern in line :

				if writing :
					# there is already a file opened, then close it
					g.close()

				g = open( path_outputs + '/trajectory_%06d' % i + '.txt' , 'a' )

				g.write( line ) 
					
				writing = True
				i = i + 1

			elif writing :
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

def intervals( xx , yy , E_min = 5 ) : 

	# compute the number of observed (O) and expected (E) counts in bins defined by
	#						np.histogram_bin_edges
	# These bins span the ecentricity range covered by the observed eccentricities y.
	# Bins are constrained to have at least 5 counts each for the chisq accuracy. 
	# Neighbour bins with small counts are merged together. It can happen that the
	# number of bins n is smaller than the recommended number of 4. However, there are 
	# no constrains on E and O, which technically allows n = 1. Obviously, the chisq 
	# accuracy will be extremely poor. Acquire movies with shorter exposure times to 
	# counter that.

	if len( yy ) != len( xx ) :
		raise AttributeError( 'intervals: x and y number of not nan values differs' )

	# use hist auto bins to find the intervals I starting from the Observed yy
	I = np.histogram_bin_edges( xx + yy , bins = 'auto' )
	#I = np.histogram_bin_edges( xx , bins = 'auto' )

	# count the Expected and Observed counts in each interval
	E = []
	O = []
	for i in range( len( I ) - 1 ) :
		
		if i == 0 :

			E.append( len( [ x for x in xx if ( ( x >= I[ i ] ) & ( x <= I[ i + 1 ] ) ) ] ) )
			O.append( len( [ y for y in yy if ( ( y >= I[ i ] ) & ( y <= I[ i + 1 ] ) ) ] ) )

		elif i == len( I ) - 1 :

			E.append( len( [ x for x in xx if ( ( x > I[ i ] ) & ( x <= I[ i + 1 ] ) ) ] ) )
			O.append( len( [ y for y in yy if ( ( y > I[ i ] ) & ( y <= I[ i + 1 ] ) ) ] ) )

		else : 

			E.append( len( [ x for x in xx if ( ( x > I[ i ] ) & ( x <= I[ i + 1 ] ) ) ] ) )
			O.append( len( [ y for y in yy if ( ( y > I[ i ] ) & ( y <= I[ i + 1 ] ) ) ] ) )

	# if there are Expected counts smaller than Emin group those interval with the neighbour interval.
	while ( ( min( E ) < E_min ) & ( len( E ) > 1 ) ) :		# it can happen that there are no values in
															# the bins identified by the observed values 
		EE = []												# yy, if so, all E will be 0 and will be 
		OO = []												# groued into one E = [ 0 ] bin. At this 
		II = [ I[ 0 ] ]										# point the while loop can stop, hence the 
															# condition len( E ) > 1 
		e = o = 0

		for i in range( len( E ) ) :
			
			e = e + E[ i ]
			o = o + O[ i ]
			
			if i == len( E ) - 1 :
				
				if ( ( e < E_min ) & ( len( EE ) > 0 ) ) :	# if all e were smaller than E_min, no value
															# has been appended to EE, hence EE[-1] does 
					EE[-1] = EE[-1] + e						# not exist. In this case, skip this option
					OO[-1] = OO[-1] + o						# append what is in e to EE (next else), and
					II[-1] = I[ i + 1 ]						# exit the while loop ( len( E ) == 1 )

				else :
					
					EE.append( e )
					OO.append( o )
					II.append( I[ i + 1 ] )
			
			elif e < E_min :
				
				continue
			
			else :
				
				EE.append( e )
				OO.append( o )
				II.append( I[ i + 1 ] )
				e = o = 0

		E = EE
		O = OO
		I = II
		
	n = len( E )

	if n == 1 : 

		print( "Warning from intervals: only one interval with <= " + str( E_min ) + " expected observations. Consider inputing trajectories with more datapoints" )

	return O , E , n , I

def chi2test( O , E , n , v ) :
	"""
	chisqu( O , E , n , v = 0 ) coputes the chi squared statistics on the number of 
	(O)bserved VS (E)xpected counts. n is the number of bins and (v) is the number of 
	constrains. 
	"""
	
	l = len( O )
	if len( E ) != l :
		raise AttributeError( 'chisq: O and E numbers differ' )
	
	df = n - v 
	try :
	
		chi  = sum( [ ( O[ i ] - E[ i ] ) ** 2 / E[ i ] for i in range( l ) ] )
		p = 1 - chi2.cdf( chi , df )
	
	except :
		
		chi = np.nan
		p = 0 

	return p , chi , df

def eccStats( t , rt , v = 1 , fimax = True , plot = False , verbose = True ) :

	filename = t.annotations()[ 'file' ]

	# compare if the eccentricity values compute in the region where the spot
	# is quantified and in the surrounding region are falling on a line close to the 
	# diagonal: y = m0 * x + c0. The H0 is that y = m0 * x + c0 is sufficient to describe 
	# the correlation between the eccentricities (i.e. the spot needs to be kept).
	if fimax :
		x , _ = ecc( t.fimax() )
		y , _ = ecc( rt.fimax() )
	else :
		x , _ = ecc( t )
		y , _ = ecc( rt )

	# remove nan if any and perform a log transform on the distribution of observed and expected values
	l = len( x )
	xx = [ np.log( x[ i ] ) for i in range( l ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] # & ( x[ i ] <= np.nanmedian( x ) ) ) ]
	yy = [ np.log( y[ i ] ) for i in range( l ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] # & ( x[ i ] <= np.nanmedian( x ) ) ) ]
	
	#xx = [ x[ i ] for i in range( l ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] # & ( x[ i ] <= np.nanmedian( x ) ) ) ]
	#yy = [ y[ i ] for i in range( l ) if ( ( x[ i ] == x[ i ] ) & ( y[ i ] == y[ i ] ) ) ] # & ( x[ i ] <= np.nanmedian( x ) ) ) ]

	O , E , n , i = intervals( xx , yy )

	if v == 1 : # constrain the sum( E ) == sum( O )
		
		if ( ( sum( E ) > 0 ) & ( sum( O ) > 0 ) ):
		
			E_constrained = [ E[ i ] * sum( O ) / sum( E ) for i in range( len( E ) ) ]
			E = E_constrained

		# else keept the E and O unchanged, no Expected or Observed values in the binning range cannot be changed.
		
	p , _ , _ = chi2test( O , E , n , v = v )
	
	# compute the regression coefficient
	_ , _ , r , _ , _ = linregress( xx , yy )


	if verbose :
	
		print( '-- ' + filename + ' --' )
		print( 'Observed counts (O): ' + str( O ) )
		print( 'Expected counts (E): ' + str( E ) )
		print( 'Intervals (I) : ' + str( i ) )
		print( 'Constrains (v) : ' + str( v ) )
		print( '-> p-value: ' + str( p ) )
		print( '-> R: ' + str( r ) )

	if plot :
		
		if not os.path.exists( 'Plots' ) :
			os.makedirs( 'Plots' )
			
		plt.figure( figsize = ( 12 , 5 ) )
	
		plt.subplot( 121 )
		plt.title( filename + '; ' + r'$R=$' + str( round( r , 3 ) ) ) 
		plt.plot( xx , yy , marker = 'o' , ls = '' )
		plt.plot( xx , xx , ls = '-' )
		plt.xlabel( r'$\log(\epsilon)$' )
		plt.ylabel( r'$\log(\epsilon_r)$' )

		plt.subplot( 122 )
		plt.title( r'$p(\chi^2 > \chi^2_0)=$' + str( round( p , 3 ) ) + '; ' + r'$n=$' + str( n ) + '; ' + r'$df=$' + str( n - v ) )
		plt.hist( xx , label = 'Expected (E): ' + r'$\log(\epsilon)$' , alpha = 0.5 , color = '#00ff00' )
		plt.hist( yy ,  label = 'Observed (O): ' + r'$\log(\epsilon_r)$' , alpha = 0.5 , color = '#ff0000' )
		plt.bar( i[1:] , E , width = [ i[j-1] - i[j] for j in range( 1 , len( i ) ) ] , align = 'edge' , color = 'none' , edgecolor = '#00ff00' ) 
		plt.bar( i[1:] , O , width = [ i[j-1] - i[j] for j in range( 1 , len( i ) ) ] , align = 'edge' , color = 'none' , edgecolor = '#ff0000' , ls = '--' ) 
		plt.legend()
		plt.xlabel( r'$\log(\epsilon)$' )
		plt.ylabel( 'Density' )
		
		plt.savefig( 'Plots/' + filename[:-4] + '.pdf' )
		plt.close()

	if p != p : # convert nan pvalues to 0 pvalues
		
		p = 0
	
	return p , r 

def ichose( tt , rtt , image_shape, image_len, pval = 0.05 , r_min = 0.5 , d0 = 10 , f0 = 0 , plot = True , fimax = True , verbose = True , v = 1 , output_rtt = False ) :
	"""
	ichose( tt , rtt , image_shape, image_len, maxit = 100 , d0 = 10 , f0 = 0 )
	select the trajectories in the trajectory list 'tt' whose spots have eccentricities that do not change if measured in
	a larger region around the spot. These eccentricity values are inputed through the trajectory list 'rtt'. 
	- image_shape, define the shape of the image so that d0 (see below) can be used
	- image_len, define the length of the image so that f0 (see below) can be used
	more stringent the spot selection
	- maxit, max number of iterations in rlm (R)
	- d0, the region on the image border where trajectories are rejected by default (too close to the image edges)
	- f0, trajectories starting before frame f0 and ending after frame image_len - f0 - 1 are rejected
	"""

	# define the list of selected trajectories that will be outputed
	output_good_tt = [] # selected trajectory list
	output_good_rtt = [] # selected trajectory list with excentricities in larger r
	output_bad_tt = [] # rejected trajectory list
	output_bad_rtt = [] # rejected trajectory list with excentricities in larger r

	l = len( tt ) 

	if not l == len( rtt ) :

		raise TypeError( 'the list of trajectories tt and rtt have different lengths' )

	for i in range( l ) :

		p , r = eccStats( tt[ i ] , rtt[ i ] , plot = plot , fimax = fimax , verbose = verbose , v = v )

		if ( ( p > pval ) & ( r > r_min ) ) : 

			mc = mean_centroid( tt[ i ] )

			if ( ( mc[ 0 ] > d0 ) & ( mc[ 1 ] > d0 ) ) :
				
				if ( ( mc[ 0 ] < ( image_shape[ 0 ] -  d0 ) ) & ( mc[ 1 ] < ( image_shape[ 1 ] - d0 ) ) ) :

					if ( ( tt[ i ].frames( 0 ) > f0 ) & ( tt[ i ].frames( len( tt[ i ] ) - 1 ) < image_len - 1 - f0 ) ) :
					
						# annotations useful to check selection parameters
						tt[ i ].annotations( 'eccentricity_pval' , str( p ) )
						tt[ i ].annotations( 'eccentricity_R' , str( r ) )

						# append good trajectory
						output_good_tt.append( tt[ i ] )
						output_good_rtt.append( rtt[ i ] )
		
		else :

			# append bad trajectory
			output_bad_tt.append( tt[ i ] )
			output_bad_rtt.append( rtt[ i ] )

	if output_rtt :
		
		return output_good_tt , output_good_rtt , output_bad_tt , output_bad_rtt

	else : # default
		
		return output_good_tt , output_bad_tt


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

