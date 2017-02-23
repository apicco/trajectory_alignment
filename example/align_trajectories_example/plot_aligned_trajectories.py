from trajalign.traj import Traj
from numpy import transpose, concatenate
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon

#function to plot the average trajectories and the 95% confidence interval
def myplot( obj , t , what , label , col , scale = 0.5 ) :
	
	x = getattr( t , '_' + what )
	x_err = getattr( t , '_' + what + '_err' )

	if x.ndim > 1 : 
		#then the attribute has more than one dimention and we are interested
		#only in the first one.
		x = x[ 0 ]
		x_err = x_err[ 0 ]

	lower_error_boundary =  transpose( 
			[ t.t() , x - 1.96 * x_err ]
			)
	upper_error_boundary =  transpose( 
			[ t.t() , x + 1.96 * x_err ] 
			)
	error_boundary = concatenate( ( 
		lower_error_boundary , upper_error_boundary[ ::-1 ] 
		) )

	error_area = Polygon( error_boundary , True , color = col , alpha = 0.3 )
	obj.add_patch( error_area )

	#plot the trajectory
	obj.plot( t.t() , x , linewidth = 1.5 , color = col , label = label )

#load the aligned trajectories
abp1 = Traj()
abp1.load( 'abp1.txt' )

sla1 = Traj()
sla1.load( 'sla1_aligned.txt' )

rvs167 = Traj()
rvs167.load( 'rvs167_aligned.txt' )

#normalise the fluorescence intensities
abp1.norm_f()
sla1.norm_f()
rvs167.norm_f()

#set trajectories start from time = 0 s

t_0 = min( concatenate(( sla1.t() , abp1.t() , rvs167.t() )) )
abp1.input_values( 't' , abp1.t() - t_0 )
sla1.input_values( 't' , sla1.t() - t_0 )
rvs167.input_values( 't' , rvs167.t() - t_0 )

#plot
f, ( trj , fi ) = plt.subplots( 2 , 1 , gridspec_kw = { 'height_ratios' : [ 2 , 1 ] } , figsize = ( 8 , 11 ) , sharex = True )

#abp1

myplot( trj , abp1 , what = 'coord' , col = '#D7110E' , label = 'Abp1' )
myplot( trj , sla1 , what = 'coord' , col = '#336CFF' , label = 'Sla1' )
myplot( trj , rvs167 , what = 'coord' , col = '#006400' , label = 'Rvs167' )

myplot( fi , abp1 , what = 'f' , col = '#D7110E' , label = 'Abp1' )
myplot( fi , sla1 , what = 'f' , col = '#336CFF' , label = 'Sla1' )
myplot( fi , rvs167 , what = 'f' , col = '#006400' , label = 'Rvs167' )

plt.subplot( trj )
plt.ylabel( 'Inward movement (' + abp1.annotations( 'coord_unit' ) + ')' , fontsize = 24 )
plt.legend( loc = 'best' )

plt.subplot( fi )
plt.ylabel( 'FI (a.u.)' , fontsize = 24 )
plt.xlabel( 'Time (' + abp1.annotations( 't_unit' ) + ')' , fontsize = 24 )
plt.legend( loc = 'best' )

f.tight_layout()
f.savefig( 'plot_aligned_trajectories.png' )

