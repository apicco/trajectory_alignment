# All the software here is distributed under the terms of the GNU General Public License Version 3, June 2007. 
# Trajalign is a free software and comes with ABSOLUTELY NO WARRANTY.
# 
# You are welcome to redistribute the software. However, we appreciate is use of such software would result in citations of 
# Picco, A., Kaksonen, M., _Precise tracking of the dynamics of multiple proteins in endocytic events_,  Methods in Cell Biology, Vol. 139, pages 51-68 (2017)
# http://www.sciencedirect.com/science/article/pii/S0091679X16301546
# 
# Author: Andrea Picco (https://github.com/apicco)
# Year: 2017

from numpy import array
from numpy import transpose
from numpy import matrix
from numpy import square
from numpy import sqrt
from numpy import sin
from numpy import sign
from numpy import cos
from numpy import insert
from numpy import NaN
from numpy import insert
from numpy import nanmean
from numpy import nanstd
from numpy import nanmax
from numpy import nanmin
from numpy import float64 
from numpy import convolve
from numpy import isclose
from numpy import isnan
from numpy import polyfit
from numpy import inf
from numpy import round
import copy as cp
import os
import re

class Traj:
	"""Trajectory OBJECT:
		traj(**annotations) -> creates a new empty trajectory. **annotations are
		an optional dictionary with entries about the trajectory (name, date,
		experiment,..). 
		A trajectory contains several attributes:
			.frames() -> array of frame indexes in the original image, if any.
			.t() -> arry of time values associated to each frame.
			.coord() -> spatial coordinate matrix of the fluorescent spot centroid
			.f() -> array of fluorescence intensities measured from the spots.
			.mol() -> array of number of molecules, if known.
			.n() -> if traj() is an average trajectory, then 
				traj().n output an array containing the number 
				of trajectories used to compute the average.
			.m2() -> is the second moment of brightness of the fluresencent 
				patches.

		.t(), .coord(), .f(),  .mol(), and .m2() have related attributes containing the 
		measured errors, if known. The errors are called with .t_err(), 
		.coord_err(), .f_err(), .mol_err(), and .m2_err() respectively
		
		PROPERTIES:

		Each attribute of the trajectory, with the exception of .name, allows 
		integer arguments to indexing elements in its array, which are then 
		returned.
		
		If .frames() and .t() are both non-empty arrays, they  must have the 
		same length. All trajectory elements that have a length different than 
		zero must have the same length than .frames() and/or .t().

		Only the non-emty attributes of a trajectory are printed as a table. The
		table header contains the names of the attributes that are printed. The 
		table is complemented with the annotations, if any.

		MODULES:
		
		.input_values(attribute_names,array,unit='') inputs array in the trajectroy 
		attibute attribute_names. If known, allows to add the unit associated to the
		attribute values.

		.start() returns the start of the trajectory if the time attribute is 
		defined.
		
		.end() returns the end of the trajectory if the time attribute is defined.

		.fill() fills attributes of missing frames with NaN 
		.frames() and .t() accordingly.

		.save(filename) saves the trajectory as txt to the filename.
		
		.rotate(angle) rotates the coordinates by 'angle' expressed in radiants.

		.translate(v): translates the coordinates of the trajectory by a vector
		v: v[0] shifts .x[0,] while v[1] shifts .x[1,].

		.load(filename,sep=None,comment_char='#',**attribute_names): loads data from a txt table.
		Data must be ordered in columns. Columns can be separated by spaces or
		tabs or comas (for .csv files). 
		The separator can be entered in sep as a string. The default for sep is None, 
		which will recognise colums separated by an arbitraty number of spaces.
		The strings starting with the commment character(s) in the comment_char
		string will be disregarded.
		Attribute_names associate the attribute and the number of the column 
		containing its data. Attribute_names can be only: "frames", "t", "x", 
		"f", "n", "mol", "m2", "t_err", "x_err", "f_err", "mol_err", "m2_err".

		example of usage of the load function:
		t = Traj() #empty trajectory
		t.load('filename.txt',frames=0,t=1,x=(2,3))

		Note that 'coord' requires two valuse and the column indexing starts 
		from 0.

		.annotations(annotation=None,string=''): output the dictionary of the annotations associate to the
		trajectory (equivalent to .__dict__()). If an annotation is inputed it changes
		the value of the annotation with string. If the annotation is not defined it defines it.
		
		EXAMPLES:

		import numpy as np

		t = Traj(name='My favorite protein', date='25-12-0001', experiment='The meaning of life') #creates a trajectory object called t
		
		t.input_values('frames',[0,1,2,3,4,5,6,7,8,9]) #input the frame indexes
		t.input_values('coord',[[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] \
				[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]]) #input the coords 
		
		print(t.frames()) #print the frame indexes
		print(t.frames(0,len(t)-1)) #print the first and last frame indexs.
		print(t.coord(0,len(t)-1)) #print the first and last coordinates.
		print(t.__dict__()) #print the annotations to the trajectory
		print(t) #print the trajectory content

		t.rotate(np.pi/2) #rotate the trajectory by pi/2 
		print(t)

		"""
	
	__slots__ = ['_annotations','_frames','_t','_coord','_f','_mol','_n','_m2', '_m3' , '_m4' , '_m5' , '_u02' , '_u20' , '_u11' , '_ecc' , '_t_err','_coord_err','_f_err','_mol_err' , '_m2_err' , '_m3_err', '_m4_err', '_m5_err', '_u02_err', '_u20_err', '_u11_err' ]
	

	def __init__(self,**annotations):

		#Trajectory main attributes 
		self._annotations = annotations
		self._frames = array([],dtype='int64')
		self._t = array([],dtype='float64')
		self._coord = array([array([]),array([])],dtype='float64') #equivalent to matrix in python3.5 provided I use @ operator
		self._f = array([],dtype='float64')
		self._mol = array([],dtype='float64')
		self._n = array([],dtype='float64')
		self._m2 = array([],dtype='float64') # second moment of brightness (f is the first)
		self._m3 = array([],dtype='float64') # third moment of brightness 
		self._m4 = array([],dtype='float64') # fourth moment of brightness 
		self._m5 = array([],dtype='float64') # fifth moment of brightness 
		self._u02 = array([],dtype='float64') # second moment of brightness along y 
		self._u20 = array([],dtype='float64') #	second moment of brightness along x
		self._u11 = array([],dtype='float64') # covariance brightness 
		self._ecc = array([],dtype='float64') # trajectory eccentricity

		#Trajectory error attributes
		self._t_err = array([],dtype='float64')
		self._coord_err = array([array([]),array([])],dtype='float64') #equivalent to matrix in python3.5 provided I use @ operator
		self._f_err = array([],dtype='float64')
		self._mol_err = array([],dtype='float64')
		self._m2_err = array([],dtype='float64')
		self._m3_err = array([],dtype='float64')
		self._m4_err = array([],dtype='float64')
		self._m5_err = array([],dtype='float64')
		self._u02_err = array([],dtype='float64')
		self._u20_err = array([],dtype='float64')
		self._u11_err = array([],dtype='float64')


	def __dict__(self):
		return self._annotations

	def __len__(self): #the number of timepoints in the trajectory
		if (len(self._t) > 0): return len(self._t)
		else: return len(self._frames)

	def __repr__( self , n0 = 0 , n1 = NaN ):
		if (len(self)) == 0 :
			output = 'The trajectory is empty!\n'
			if (len(self._annotations)):
				output += '#' + '-' * len(output) + '\n'
				for name, item in self._annotations.items():
					output += '# ' + name + ': ' + str( item )+ '\n'
			return output
		else :
			table = []
			names = []
			output = '#'
			for s in self.__slots__[1:]:
				x = getattr(self,s)
				# old, deprecated because not all attributes are to be represented in a table <- if (x.shape[x.ndim-1] > 0):
				if (x.shape[x.ndim-1] == len( self ) ):
					if s in ('_coord','_coord_err'):
						#x coord
						table.append(x[0])

						try :
							names.append('x' + s[6:] + ' (' + self._annotations['coord_unit'] + ')')
						except :
							names.append('x' + s[6:])

						#y coord	
						table.append(x[1])
						
						try :
							names.append('y' + s[6:] + ' (' + self._annotations['coord_unit'] + ')')
						except:
							names.append('y' + s[6:])
					else: 
						table.append(x)
# old						if s[1:] + '_unit' in self._annotations.keys():
# old							if len(self._annotations[ s[1:] + '_unit' ]) > 0:
# old								names.append( s[1:] + ' (' + self._annotations[s[1:] + '_unit' ] + ')' )
# old							else: 
# old								names.append(s[1:])
# old						else:
# old							names.append(s[1:])
						try :
							names.append( s[1:] + ' (' + self._annotations[s[1:] + '_unit' ] + ')' )
						except :
							names.append(s[1:])
			#find the best column width for the table
			table_col_width = max(len(str(elmnt)) \
					for row in table for elmnt in row) + 2
			name_width = max(len(elmnt) for elmnt in names) + 2
			col_width = max(table_col_width,name_width)
			#print header
			output += str(names[0]).rjust(col_width-1)#the line begins with a '#'
			for name in names[1:]:
				output += str(name).rjust(col_width)
			output += '\n'
			#print table
			row_ID = 0
			if not n1 == n1 : n1 = len( self )
			for row in transpose(table):
				if ( row_ID >= n0 ) and ( row_ID < n1 ) :
					output_row = ''
					for element in row:
						output_row += str(element).rjust(col_width)
					output += output_row+'\n'
				row_ID += 1
			#print annotations
			output += '#'+'-'*(len(names)*col_width-1)+'\n'#nice separator
			for name, item in self._annotations.items():
				output += '# '+name+': '+ str( item )+'\n'
			return output


	#Getters
	def slots( self ) :

		return self.__slots__[1:]

	def frames(self,*items):
		if (len(items)==0): return self._frames
		elif len( items ) == 1 : return self._frames[ items ]
		else: 
			try:
				return(self._frames[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().frames are out of bounds')

	def t(self,*items):
		if (len(items)==0): return self._t
		elif len( items ) == 1 : return self._t[ items ]
		else: 
			try:
				return(self._t[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().t are out of bounds')
	
	def coord(self,*items):
		new_items = [] #items are enters as a tuple, and should be converted as list
		for i in items:
			if isinstance(i, int):
				new_items.append(i)
			else:
				while isinstance(i,list):
					if len(i) == 1:
						i = i[0] #remove nested lists of lists of length 1
					else:
						i = tuple(i) #exit the while loop
				if isinstance(i,int) :
					new_items.append(i)
				else :
					for k in i:
						new_items.append(k)
		if (len(new_items)==0): return self._coord
		else: 
			try:
				output=self._coord[:,[item for item in new_items]]
				return(output)
			except IndexError:
				print('Indexes in Traj().coord are out of bounds')

	def f(self,*items):
		if (len(items)==0): return self._f
		elif len(items) == 1 : return self._f[ items ]
		else: 
			try:
				return(self._f[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().f are out of bounds')

	def mol(self,*items):
		if (len(items)==0): return self._mol
		elif len(items) == 1 : return self._mol[ items ]
		else: 
			try:
				return(self._mol[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().mol are out of bounds')

	def n(self,*items):
		if (len(items)==0): return self._n
		elif len(items) == 1 : return self._n[ items ]
		else: 
			try:
				return(self._n[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().n are out of bounds')

	def m2(self,*items):
		if (len(items)==0): return self._m2
		elif len(items) == 1 : return self._m2[ items ]
		else: 
			try:
				return(self._m2[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m2 are out of bounds')

	def m3(self,*items):
		if (len(items)==0): return self._m3
		elif len(items) == 1 : return self._m3[ items ]
		else: 
			try:
				return(self._m3[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m3 are out of bounds')

	def m4(self,*items):
		if (len(items)==0): return self._m4
		elif len(items) == 1 : return self._m4[ items ]
		else: 
			try:
				return(self._m4[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m4 are out of bounds')

	def m5(self,*items):
		if (len(items)==0): return self._m5
		elif len(items) == 1 : return self._m5[ items ]
		else: 
			try:
				return(self._m5[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m5 are out of bounds')

	def u02(self,*items):
		if (len(items)==0): return self._u02
		elif len(items) == 1 : return self._u02[ items ]
		else: 
			try:
				return(self._u02[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u02 are out of bounds')

	def u20(self,*items):
		if (len(items)==0): return self._u20
		elif len(items) == 1 : return self._u20[ items ]
		else: 
			try:
				return(self._u20[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u20 are out of bounds')

	def u11(self,*items):
		if (len(items)==0): return self._u11
		elif len(items) == 1 : return self._u11[ items ]
		else: 
			try:
				return(self._u11[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u11 are out of bounds')

	def t_err(self,*items):
		if (len(items)==0): return self._t_err
		elif len(items) == 1 : return self._t_err[ items ]
		else: 
			try:
				return(self._t_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().t_err are out of bounds')

	def coord_err(self,*items):
		new_items = []
		for i in items:
			if isinstance(i, int):
				new_items.append(i)
			else:
				while isinstance(i,list):
					if len(i) == 1:
						i = i[0] #remove nested lists of lists of length 1
					else:
						i = tuple(i) #exit the while loop
				if isinstance(i,int) :
					new_items.append(i)
				else :
					for k in i:
						new_items.append(k)
		if (len(new_items)==0): return self._coord_err
		else: 
			try:
				output=self._coord_err[:,[item for item in new_items]]
				return(output)
			except IndexError:
				print('Indexes in Traj().coord_err are out of bounds')

	def f_err(self,*items):
		if (len(items)==0): return self._f_err
		elif len(items) == 1 : return self._f_err[ items ]
		else:
			try:
				return(self._f_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().f_err are out of bounds')
	
	def mol_err(self,*items):
		if (len(items)==0): return self._mol_err
		elif len(items) == 1 : return self._mol_err[ items ]
		else: 
			try: 
				return(self._mol_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().mol_err are out of bounds')

	def m2_err(self,*items):
		if (len(items)==0): return self._m2_err
		elif len(items) == 1 : return self._m2_err[ items ]
		else:
			try:
				return(self._m2_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().f_err are out of bounds')
	
	def m3_err(self,*items):
		if (len(items)==0): return self._m3_err
		elif len(items) == 1 : return self._m3_err[ items ]
		else: 
			try:
				return(self._m3_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m3_err are out of bounds')

	def m4_err(self,*items):
		if (len(items)==0): return self._m4_err
		elif len(items) == 1 : return self._m4_err[ items ]
		else: 
			try:
				return(self._m4_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m4_err are out of bounds')

	def m5_err(self,*items):
		if (len(items)==0): return self._m5_err
		elif len(items) == 1 : return self._m5_err[ items ]
		else: 
			try:
				return(self._m5_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().m5_err are out of bounds')

	def u02_err(self,*items):
		if (len(items)==0): return self._u02_err
		elif len(items) == 1 : return self._u02_err[ items ]
		else: 
			try:
				return(self._u02_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u02_err are out of bounds')

	def u20_err(self,*items):
		if (len(items)==0): return self._u20_err
		elif len(items) == 1 : return self._u20_err[ items ]
		else: 
			try:
				return(self._u20_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u20_err are out of bounds')

	def u11_err(self,*items):
		if (len(items)==0): return self._u11_err
		elif len(items) == 1 : return self._u11_err[ items ]
		else: 
			try:
				return(self._u11_err[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().u11_err are out of bounds')

	def ecc(self,*items):
		if (len(items)==0): return self._ecc
		elif len(items) == 1 : return self._ecc[ items ]
		else: 
			try:
				return(self._ecc[[item for item in items]])
			except IndexError:
				print('Indexes in Traj().ecc are out of bounds')


	def extract(self,*items):
		
		"""
		extract(*items): extracts the rows of the trajectory defined by *items and happend annotations. Some examples:
		.extract(10,11,15) extracts rows 10, 11 and 15.
		.extract([10,11,15]) extracts rows 10, 11 and 15.
		.extract(range(10,20)) extracts all the rows from 10 to 19. 
		"""

		if (len(items)==0): 
			raise IndexError('Please, specify the values you want to extract from the trajectory')
		else:
			new_items = [] #items are enters as a tuple, and should be converted as list
			for i in items:
				if isinstance(i, int):
					new_items.append(i)
				else:
					while isinstance(i,list):
						if len(i) == 1:
							i = i[0] #remove nested lists of lists, which have length 1
						else:
							i = tuple(i) #exit the while loop
					if isinstance(i,int) :
						new_items.append(i)
					else:
						for k in i:
							new_items.append(k)
			#create the output trajectory
			if 'range' in self.annotations().keys():
				output = Traj(range = self.annotations()['range']+' then '+str(new_items))
			else:
				output = Traj(range = str(new_items))
			#inherit the annotations
			for a in self.annotations().keys():
				if a != 'range':
					output.annotations(a,self._annotations[a])
			try:
				for a in self.attributes():
					if a == 'frames' : 
						output.input_values(a,self.frames(new_items))
					if (a == 't') : 
						output.input_values(a,self.t(new_items))
					if a == 'coord' : 
						output.input_values(a,self.coord(new_items))
					if a == 'f' : 
						output.input_values(a,self.f(new_items))
					if a == 'mol' : 
						output.input_values(a,self.mol(new_items))
					if a == 'n' : 
						output.input_values(a,self.n(new_items))
					if a == 'm2' : 
						output.input_values(a,self.m2(new_items))
					if a == 'ecc' : 
						output.input_values(a,self.ecc(new_items))
					if a == 't_err' : 
						output.input_values(a,self.t_err(new_items))
					if a == 'coord_err' : 
						output.input_values(a,self.coord_err(new_items))
					if a == 'f_err' : 
						output.input_values(a,self.f_err(new_items))
					if a == 'mol_err' : 
						output.input_values(a,self.mol_err(new_items)) 
					if a == 'm2_err' : 
						output.input_values(a,self.m2_err(new_items)) 
			except IndexError:
				print('Indexes in range are out of bounds')

		return output

	def head( self , n = 10 ) :

		"""
		.head( n = 10 ) returns the first n rows in the trajectory 
		"""

		print( self.__repr__( n1 = n ) )

	def tail( self , n = 10 ) :
		
		"""
		.tail( n = 10 ) returns the last n rows in the trajectory 
		"""

		print( self.__repr__( n0 = len( self ) - n ) )

	def fimax( self , filter = [ 1 ] ):

		"""
		.fimax( self , filter = [ 1 ] ) : extracts a trajectory that stops at the max of the fluorescence intensity, included.
		'filter' defines the  filter used to smooth the fluorescence intensity profile. Default is no filter ( filter = [ 1 ] ).
		"""

		output = cp.deepcopy( self )
		
		#check that the trajectory has a fluorescence intensity attribute which is not empty
		if not len( self.f() ) :

			raise AttributeError('The fluorescence intensity attribute is empty!')

		if not len( self.t() ) :
			
			raise AttributeError('The time attribute is empty!')

		#running mean over fi. The same mean is run over the row numbers to find the 
		#row at which the max in fi is.
		not_nan = [ i for i in range( len( self ) ) if self.f( i ) == self.f( i ) ]
		t_to_convolve = self.extract( not_nan )
		
		fi = array( convolve( t_to_convolve.f() , filter , 'valid' ) )
		times = array( convolve( t_to_convolve.t() , filter , 'valid' ) )

		time_where_fi_max_is = times[ min( fi.argmax() + 1 , len( fi ) - 1 ) ]

		output.end( time_where_fi_max_is )

		#create annotation
		output.annotations( 'fimax' , 'TRUE' )

		return output

	def msd( self , scale = 1 ) :

		#check that the attribute .coord is not empty
		if len( self.coord() ) != 2 : 

			raise AttributeError( 'The size of the .coord() attribute does not match expectation, does your trajectory have a x and y coordinates?' ) 
		
		# check that the annotation delta_t is not emplty
		if not self.annotations()[ 'delta_t' ] :
			
			raise AttributeError( 'msd: The trajectory has no annotation delta_t!' ) 

		else :

			# fill the trajectory, it is important that missing gaps are filled
			self.fill()
		
			# msd units:
			if ( ( scale == 1 ) & ( self.annotations()[ 'coord_unit' ] == 'pxl' ) ):
	
				raise AttributeError( 'coordinates are in pxl units, add a scaling factor' )
			
			elif ( ( scale != 1 ) & ( self.annotations()[ 'coord_unit' ] != 'pxl' ) ) :
			
				raise AttributeError( 'coordinates are in units, you do not need a scaling factor' )
		
			m = [ 0 ]
			sem = [ inf ]
		
			l = len( self.coord()[ 0 ] )
			ss = 1  #initiate the step size

			while( ss < l ) :

				x1 = self.coord()[ 0 ][ ss : ] * scale
				x2 = self.coord()[ 0 ][ : l - ss ] * scale
				
				y1 = self.coord()[ 1 ][ ss : ] * scale
				y2 = self.coord()[ 1 ][ : l - ss ] * scale

				d = [ ( x1[ i ] - x2[ i ] ) ** 2  + ( y1[ i ] - y2[ i ] ) ** 2 for i in range( l - ss ) ] 

				m.append( nanmean( d ) )
				sem.append( nanstd( d ) / sqrt( sum( ~isnan( d ) ) ) ) #sem, dividing by the square root of the number of not nan items

				ss = ss + 1 
			
			return array( [ 
				[ i * float( self.annotations()[ 'delta_t' ] ) for i in range( len( m ) ) ] ,
				m , 
				sem ] 
				)

	def msdfit( self , sel = None , scale = 1 , deg = 2 , return_data = False ) :

		if sel == None :
			sel = len( self )

		# output dt , msd , msd_err into m
		m = self.msd( scale ) 
		t = m[ 0 ][ 0 : sel ]
		y = m[ 1 ][ 0 : sel ]
		y_err = m[ 2 ][ 0 : sel ]
		
		# compute the msd polyfit of deg 
		try: 

			p , cov = polyfit( t , y , w = 1/y_err , deg = deg , cov = True )
		
		except :

			p , cov = polyfit( t , y , deg = deg , cov = True )

		v = [ sqrt( p[ 0 ] ) , sqrt( cov[ 0 , 0 ] ) / ( 2 * sqrt( p[ 0 ] ) ) ]

		D = [ p[ 1 ] / 4 , cov[ 1 , 1 ] / 4 ]
	
		if return_data :
	
			data = { 'x' : t , 'y' : y , 'err' : y_err , 'p' : p , 'cov' : cov }
			
			return v , D , data

		else : 
			
			return v , D

	#Setters
	#Input values in the Traj object as arrays. Array length must be equal to the length of frames and time
	def input_values(self,name,x,unit=''):
		"""
		input_values(attribute_names,array,unit='') inputs array in the trajectroy 
		attibute attribute_names. If known, allows to add the unit associated to the
		attribute values.
		"""
		if ('_'+name in self.__slots__[1:]):

			if ((name=='frames') & (len(self._t)==0)):
				#copute the extent of gaps between frames (in general it is >= 1). If a 
				#gap is negative it means that the chronological order of the frames
				#is erroneous.
				frame_gaps = [ x[ i + 1 ] - x[ i ] for i in range( 0 , len( x ) - 1 ) ]
			
				if len( x ) == 1 : #if there is only one frame it is not possible to compute frame_gaps
					setattr(self,"_"+name,array(x,dtype='int64')) # add the frames; no time present yet 
				elif min( frame_gaps ) > 0: 
					setattr(self,"_"+name,array(x,dtype='int64')) # add the frames; no time present yet 
				else: 
					raise AttributeError('The chronological order of the frames is wrong')

				if ( unit != '' ) : 
					raise AttributeError( 'Frames do not have units' )
				
			elif ((name=='t') & (len(self._frames)==0)):
			
				setattr(self,"_"+name,array(x,dtype='float64')) # add the time; no frames present yet
			
				if ( unit != '' ):

					self._annotations['t_unit'] = unit

			elif ( name == 'coord' ):

				if (len(x[0])==(len(self._frames) | len(self._t))):

					setattr(self,"_"+name,array(x,dtype='float64')) # add the coords
					
					if ( unit != '' ):
						self._annotations['coord_unit'] = unit
				
				else:

					raise AttributeError('The input length of the attribute does not match frames or time array of non-zero length')
			
			elif ( name == 'coord_err' ):

				if (len(x[0])==(len(self._frames) | len(self._t))):

					setattr(self,"_"+name,array(x,dtype='float64')) # add the coords

					# in the case the user wants to redefine the coord unit he needs to make sure that coord_units has
					# not yet been assigned in .coord() 
					if ( unit != '' ) :

						if ('coord_unit' in self._annotations.keys()):
							
							if ( self._annotations['coord_unit'] == '' ): 
								
								self._annotations['coord_unit'] = unit
							
							else :
	
								raise AttributeError( 'coord_unit has already been assigned for .coord(). You cannot overwrite it!' )
	
						else: self._annotations['coord_unit'] = unit
				else:

					raise AttributeError('The input length of the attribute does not match frames or time array of non-zero length, or coord_err is not defined for this trajectory')

			elif (len(x)==(len(self._frames) | len(self._t))): 

				if (name=='frames'): setattr(self,"_"+name,array(x,dtype='int64')) #add frames array checking it is as long as frames or times
				else: setattr(self,"_"+name,array(x,dtype='float64')) #add array checking it is as long as frames or times
				if len(unit) > 0:
					self._annotations[name+' unit'] = unit
			else:
				raise AttributeError('The input length of the attribute does not match frames or time array of non-zero length')
		else:
			raise AttributeError('The attribute name does not match the allowd attributes of the trajectory class. Choose one among: \'frames\',\'t\',\'t_err\',\'coord\',\'coord_err\',\'f\',\'f_err\',\'n\'')

	def norm_f(self):
		"""
		.norm_f(): normalises the fluorescence intensities .f() between 0 and 1.
		"""
	
		# f_min and f_max need to be defined at the beginning 
		# otherwise they will change for error calculation if 
		# the if loop is entered
		
		f_min = nanmin( self._f )
		f_max = nanmax( self._f )

		self._f = ( self._f - f_min ) / ( f_max - f_min )

		#if the attribute _f_err is not empty, then propagate the errors accordingly 
		if ( self._f_err.shape[ self._f_err.ndim - 1 ] > 0 ) :
			self._f_err = self._f_err / ( f_max - f_min )
			#the aim of this function is only to rescale the fluorescence intensity
			#between 0 and 1, hence we do no propagate the error of the max(self._f)
			#and min(self._f).

	def scale_f(self , v = 1 ):
		"""
		.scale_f(): scale the fluorescence intensities so that their integral is an arbitrary value 'v'. Default is 1.
		"""
		
		S = sum( self._f )
		self._f = self._f * v / S
		#if the attribute _f_err is not empty, then propagate the errors accordingly 
		if ( self._f_err.shape[ self._f_err.ndim - 1 ] > 0 ) :
			self._f_err = self._f_err * v / S 
			#the aim of this function is only to rescale the fluorescence intensity
			#so that all fluorescence intensiteis meet the same integral, 
			#hence we do no propagate the error of the S.

	def n_mol( self , N , N_err ) :
		
		"""
		.n_mol( N , N_err ) compute the number of molecules in the trajectory using the flurescence intensity profile 
		and the N number of molecules per spot. It also propagate the errors given the N_err associated to N.
		"""

		F = ( self.f() - nanmin( self.f() ) ) #/ ( nanmax( self._f ) - nanmin( self._f ) ) <- redundant constant because F is divided by nanmean( F ), which has the same constant
		M = nanmean( F )
		self.input_values( 'mol' , F * N / M )
		self.input_values( 'mol_err' , sqrt( \
				( F * N_err / M ) ** 2 +\
				( self.f_err() * N * ( M - F / len( F ) ) / M ** 2 ) ** 2 +\
				( self.f_err( self.f().tolist().index( nanmin( self.f() ) ) ) * N * ( M - F ) / M **2 ) ** 2 )
				)

	def rotate( self , angle , angle_err = 0):
		"""
		rotate(angle): rotates the coordinated of the trajectory \
				by an angle in radiants.
		"""

		R = matrix( [[ cos( angle ) , - sin( angle ) ] , [ sin( angle ) , cos( angle ) ]] , dtype = 'float64' ) 
		sR = square( R )
		
		self._coord = array( R @ self._coord )

		#if the attribute _coord_err is not empty, then propagate the errors accordingly 
		if ( self._coord_err.shape[ self._coord_err.ndim - 1 ] > 0 ) :

			self._coord_err = array( sqrt( 
					sR @ square( self._coord_err ) + \
							square( angle_err * sqrt( 1 - sR ) @ matrix([[ 1 , 0 ] , [ 0 , -1 ]] ) @ self._coord )
					))
		
		elif angle_err > 0 : #if there is not attribute _coord_err, but there is an error then
			self.insert_values( 'coord_err' , array( square( angle_err * sqrt( 1 - sR ) @ matrix([[ 1 , 0 ] , [ 0 , -1 ]] ) @ self.coord() )
					) )

	def center_mass(self):
		"""
		center_mass(): outputs the trajectory center of mass
		"""

		return( array( [ nanmean( self._coord[0,] ), nanmean( self._coord[1,] ) ] ))

	def translate( self , v , v_err = ( 0 , 0) ):
		"""
		translate(v): translates the coordinates of the trajectory \
				by a vector v: v[0] shifts .x[0,] while v[1] shifts \
				.x[1,].
		"""

		self._coord[ 0 , ] = self._coord[ 0 , ] + v[ 0 ]
		self._coord[ 1 , ] = self._coord[ 1 , ] + v[ 1 ]
		if len( v_err ) != 2 :
			raise AttributeError('The error must be a vector of length 2')
		else :
			if ( v_err[ 0 ] != 0 ) | ( v_err[ 1 ] != 0 ) :
				#if the attribute _coord_err is not empty, then propagate the errors accordingly 
				if ( self._coord_err.shape[self._coord_err.ndim-1] > 0 ) :
					self._coord_err[ 0 , ] = sqrt( self._coord_err[ 0 , ] ** 2 + v_err[ 0 ] ** 2 )
					self._coord_err[ 1 , ] = sqrt( self._coord_err[ 1 , ] ** 2 + v_err[ 1 ] ** 2 )
				else :
					setattr( self , '_coord_err' , array( [\
							 [ v_err[ 0 ] ] * ( len( self )  ),
							 [ v_err[ 1 ] ] * ( len( self )  )
								] , dtype = 'float64' ) )
	
	def lag(self,shift):
		"""
		lag(shift): shifts the time of the trajectory by 'shift', in the trajectory units. Shift is an integer that measure the number of time intervals, or frames, the trajectory has to be shifted. 
		"""
		if isinstance(shift,int):
			if len(self._t) == 0:
				raise AttributeError('There is no time to be shifted')
			elif 'delta_t' in self._annotations.keys():
				self._t += shift * float(self._annotations['delta_t'])
				return self._t
			else :
				print("Waring: lag() estimates the delta_t from the trajectory time attribute")
				delta_t = min(self._t[1:]-self._t[0:(len(self._t)-1)])
				self._t += shift * delta_t
				return self._t
		else :
			raise TypeError('shift in lag() must be integer')

	def start(self,t=None):

		"""
		start(t=None): the start time of the trajectory. If t is specified
		the trajecotry points starting from t are extracted. 
		"""
		
		if len(self._t) > 0:
			
			if 'delta_t' in self._annotations.keys():
				delta_t = float64( self._annotations['delta_t'] )
			else: 
				delta_t = min( self._t[1:] - self._t[ 0 : ( len(self._t) - 1 ) ] )
		
			if t is None:
				return self._t[0]
			elif ( 
					( ( t  > self._t[0] ) | isclose( t , self._t[0] ) | isclose( self._t[0] , t ) ) & 
					( ( t < self._t[len(self)-1] ) | isclose( t , self._t[len(self)-1] ) | isclose( self._t[len(self)-1] , t ) )
						): #check wheter t is comprised between self._t[0] and self._t[len(self)-1]. The two isclose are needed because in rare cases isclose order of argumants can lead to different results, see numpy documentation
				new_t = array([ i for i in self._t if ( (i > t) | isclose( i , t ) | isclose( t , i ) ) ] )
				new_start = self._t.tolist().index(new_t[0])
				self._t = new_t
				for attribute in self.attributes():
					if attribute == 'coord':
						self._coord = array([\
								self._coord[0][new_start:],
								self._coord[1][new_start:]
								])
					elif attribute == 'coord_err':
						self._coord_err = array([\
								self._coord_err[0][new_start:],
								self._coord_err[1][new_start:]
								])
					elif  attribute != 't':
						x = getattr(self,'_'+attribute)
						setattr(self,'_'+attribute,x[new_start:])	
			elif t > self._t[len(self)-1]:
				raise AttributeError('t is larger than the trajectory last time point')
			elif t < self._t[0]:

				def float_range(x , y , step): 

					float_range_output = [];
					while ( ( x > y ) | isclose( x , y ) | isclose( y , x ) ):
						float_range_output.append( x )
						x -= step
					return list( reversed( float_range_output ) )
				
				new_t = float_range( self._t[ 0 ] - delta_t , t , delta_t )

				self._t = insert( self._t , 0 , new_t )
				for attribute in self.attributes():
					if attribute == 'coord':
						self._coord = array([\
								insert(self._coord[0],0,[float('NaN')]*len(new_t)),
								insert(self._coord[1],0,[float('NaN')]*len(new_t))
								])
					elif attribute == 'coord_err':
						self._coord_err = array([\
								insert(self._coord_err[0],0,[float('NaN')]*len(new_t)),
								insert(self._coord_err[1],0,[float('NaN')]*len(new_t))
								])
					elif attribute == 'frames':
						new_frames = [self._frames[0] - f for f in range(len(new_t),0,-1)]
						self._frames = insert(self._frames,0,new_frames)
					elif  attribute != 't':
						x = insert(getattr(self,'_'+attribute),0,[float('NaN')]*len(new_t))
						setattr(self,'_'+attribute,x)	
			else:
				raise IndexError('The time attribute is empty')

	def end(self,t=None):
	
		"""
		end(t=None): the end time of the trajectory. If t is specified
		the trajecotry points ending before t are extracted. 
		"""

		if len(self._t) > 0:
			if 'delta_t' in self._annotations.keys():
				delta_t = float( self._annotations['delta_t'] )
			else: 
				delta_t = min(self._t[1:]-self._t[0:(len(self)-1)])
			if t is None:
				return self._t[len(self)-1]
			elif ( 
					( ( t  > self._t[0] ) | isclose( t , self._t[0] ) | isclose( self._t[0] , t ) ) & 
					( ( t < self._t[len(self)-1] ) | isclose( t , self._t[len(self)-1] ) | isclose( self._t[len(self)-1] , t ) )
						) : #in rare cases isclose order of argumants can lead to different results, see numpy documentation
				self._t = array([ i for i in self._t if ( (i < t) | isclose(i,t) ) ])
				new_end = len(self._t)
				for attribute in self.attributes():
					if attribute == 'coord':
						self._coord = array([\
								self._coord[0][0:new_end],
								self._coord[1][0:new_end]
								])
					elif attribute == 'coord_err':
						self._coord_err = array([\
								self._coord_err[0][0:new_end],
								self._coord_err[1][0:new_end]
								])
					elif  attribute != 't':
						x = getattr(self,'_'+attribute)
						setattr(self,'_'+attribute,x[0:new_end])	
			elif t < self._t[0]:
				raise AttributeError('t is smaller than the trajectory first time point')
			elif t > self._t[len(self)-1]:
				l = len(self) #store the length of the self, that will change when\
						#new time points will be added. In fact, len(self) measures\
						#the self._t length if present, which is the first attribute\
						#of the trajectory to be changed
				def float_range(x,y,step): 
					
					float_range_output = [];
					while ( ( x < y ) | isclose( x , y ) | isclose( y , x ) ):
						float_range_output.append(x)
						x += step
					return float_range_output
				new_t = float_range(self._t[l-1]+delta_t,t,delta_t)
				self._t = insert(self._t,l,new_t)
				for attribute in self.attributes():
					if attribute == 'coord':
						self._coord = array([\
								insert(self._coord[0],l,[float('NaN')]*len(new_t)),
								insert(self._coord[1],l,[float('NaN')]*len(new_t))
								])
					elif attribute == 'coord_err':
						self._coord_err = array([\
								insert(self._coord_err[0],l,[float('NaN')]*len(new_t)),
								insert(self._coord_err[1],l,[float('NaN')]*len(new_t))
								])
					elif attribute == 'frames':
						new_frames = [self._frames[l-1] + f + 1 for f in range(0,len(new_t))]
						self._frames = insert(self._frames,l,new_frames)
					elif  attribute != 't':
						x = insert(getattr(self,'_'+attribute),l,[float('NaN')]*len(new_t))
						setattr(self,'_'+attribute,x)	

		else:
			raise IndexError('The time attribute is empty')

	def lifetime(self,round=2):
		"""
		lifetime(round=2) computes the lifetime of the trajectory and rounds \
				it to 'round' number of digits.
		"""
		return self.end() - self.start()

	def time(self,delta_t,unit):
		
		"""
		time(delta_t,unit): assigns the time attribute based on the frame numbering
		given the known interval delta_t between frames
		"""
	
		if (len(self._frames) > 0 & len(self._t) == 0):
			self._t=self._frames*delta_t
			self._annotations['t_unit'] = unit
			self._annotations['delta_t'] = str(delta_t)
		else:
			if (len(self._frames) == 0) :raise AttributeError('The frames attribute was not  defined and is needed to compute the time()')
			if (len(self._t) > 0) :raise AttributeError('The time attribute is already defined')
	
	def save(self,file_name):
		if file_name[len(file_name)-3:] != 'txt' :
			file_name += '.txt'
		with open(file_name,'w') as f:
			f.write(repr(self))
		f.close()
	
	def load2( self , file_name , sep=None , coord_unit = '' , t_unit = '' , comment_char='#' , **attrs ):

		# define the columns that will be exctracted from the file based on the number of attributes requested to load
		
		columns = {}

		for a in attrs.keys() :

			if '_' + a in self.__slots__[ 1 : ] :

				if ( a == 'coord' ) | ( a == 'coord_err' ) :

					columns[ a + '_x' ] = []
					columns[ a + '_y' ] = []

				else :

					columns[ a ] = []
					
			else :

				raise TypeError( 'The attribute ' + a + ' is ill defined does not have a correspondance in self.__slots__' )

		#load the file content

		with open( file_name , 'r' ) as file :
	
			for line in file :
	
				# load Annotations
				if comment_char in line :
				
					line_elements = line.split( sep )
					
					# if the first element of the line is the comment_char, as for the Annotations
					if line_elements[ 0 ][ 0 : len( comment_char ) ] == comment_char :
		
						# and there is more than one element in the line,
						if len( line_elements ) > 1 :
						
							# and the second element is the name of an annotation, therefore with no-0 length
							if len( line_elements[ 1 ] ) > 0 :
							
								# and ending with ":"
								if line_elements[ 1 ][ -1 ] == ":" :
		
									annotation_name = line_elements[ 1 ][ : -1 ]
									
									annotation = str("")
									
									for i in range( 2 , len( line_elements ) - 1 ) :
			
										annotation = annotation + str( line_elements[ i ] ) + " " 
			
									annotation = annotation + str( line_elements[ len( line_elements ) - 1 ] ) #the last element has no sep following

									# assign Annotations
									
									self.annotations( annotation_name , annotation )

				# load Trajectory

				else :
					
					line_elements = line.split( sep )
				
					if len( line_elements ) > 0 : 

						for a in attrs.keys() :
	
							if ( a == 'coord' ) | ( a == 'coord_err' ) :
	
								columns[ a + '_x' ].append( float(  line_elements[ attrs[ a ][ 0 ] ] ) )
								columns[ a + '_y' ].append( float(  line_elements[ attrs[ a ][ 1 ] ] ) )
							
							elif ( a == 'frames' ) :
								
								columns[ a ].append( int( line_elements[ attrs[ a ] ] ) )
	
							else :
								
								try :
							
									columns[ a ].append( float(  line_elements[ attrs[ a ] ] ) )
								
								except : 
								
									raise TypeError( 'The attrs ' + a + ' expects only one value' )

		# control that the inputs are correct, namely that if t and coord are inputed then the user specifies their units
	
		if ( 'coord' in attrs.keys() ) & ( len( coord_unit ) == 0 ) :
				
			try :
				
				coord_unit = self.annotations()[ 'coord_unit' ]

			except :
				
				raise AttributeError( 'Please, specify the coordinate unit \'coord_unit\' in load' )
				
		if ( 't' in attrs.keys() ) & ( len( t_unit ) == 0 ) :
				
			try :
				
				t_unit = self.annotations()[ 't_unit' ]

			except :
				
				raise AttributeError( 'Please, specify the t unit \'t_unit\' in load' )
				

		# assign Trajectory values to Traj object

		for a in attrs.keys() :

			if ( a == 'coord' ) | ( a == 'coord_err' ) :
					
				self.input_values( a , [ columns[ a + '_x' ] ,  columns[ a + '_y' ] ] , unit = coord_unit )
			
			elif a == 't' :
					
					self.input_values( a , columns[ a ] , unit = t_unit )

			elif a == 'frames' :

				try: 
				
					self.input_values( a , columns[ a ] )
				
				except:
			
					raise AttributeError('.load_data: chronological disorder in trajectory "' + file_name + '".')

			else : 

				self.input_values( a , columns[ a ] )

	def load(self,file_name,sep=None,comment_char='#',**attrs):
		"""
		.load(file_name,sep=None,comment_char='#',**attribute_names): loads data from a txt table.
		Data must be ordered in columns. Columns can be separated by spaces or
		tabs or comas (for .csv files). 
		The separator can be entered in sep as a string. The default for sep is None, 
		which will recognise colums separated by an arbitraty number of spaces.
		The strings starting with the commment character(s) in the comment_char
		string will be disregarded.
		Attribute_names associate the attribute and the number of the column 
		containing its data and additional annotations to be added to the 
		trajectory. Attribute_names used to load the data can be only:
		"frames", "t", "x", "f", "n", "mol","t_err", "x_err", "f_err", "mol_err".

		example of usage of the load function:
		t = Traj() #empty trajectory
		t.load('filename.txt',frames=0,t=1,x=(2,3),description='A nice trajectory', 
		date = 'Yesterday')

		Note that 'coord' requires two values and the column indexing starts 
		from 0.
		"""

		output = {}
		for a in [a for a in attrs.keys() if '_'+a in self.__slots__[1:]]:
			if (a == 'coord') | (a == 'coord_err'):
				output[a] = [[],[]]
			else:
				output[a] = []
		
		# annotate the file_name
		self.annotations( 'file' , file_name )	

		with open( file_name , 'r' ) as file:
			
			for line in file:
				line_elements = line.split( sep )
		
				if len( line_elements ) > 0:
					#if the attributes are empty, the first commented line is the one with the column names
					if ( ( len( attrs.keys() ) == 0 ) & ( line_elements[ 0 ][ 0:len( comment_char ) ] == comment_char ) ) :
						attrs = {}
						i = 0
						for e in line_elements :
							if e[0] not in ('#','(','y'):
								if e[0] == 'x' :
									attrs[ 'coord'+e[1:] ] = (i,i+1) #ASSUME x y ON TWO CONSECUTIVE COLUMNS. NEED TO BE STRENGTHENED
									i += 2
								else :
									attrs[ e ] = i
									i += 1
						for a in [a for a in attrs.keys() if '_'+a in self.__slots__[1:]]:
							if (a == 'coord') | (a == 'coord_err'):
								output[a] = [[],[]]
							else:
								output[a] = []
					elif (( len(attrs.keys()) > 0 ) & ( line_elements[0][0:len(comment_char)] != comment_char )):
						for a in [a for a in attrs.keys() if '_'+a in self.__slots__[1:]]:
							try:
								if (a == 'coord') | (a == 'coord_err'):
									output[a][0].append(float(line_elements[attrs[a][0]]))
									output[a][1].append(float(line_elements[attrs[a][1]]))
								elif (a == 'frames'):
									output[a].append(int(float(line_elements[attrs[a]])))
								else :
									output[a].append(float(line_elements[attrs[a]]))
							except:
								raise TypeError('The comment_char might be ill-defined (default is "#") or the column numbering is wrong.')
					elif  line_elements[ 0 ][ 0:len( comment_char ) ] == comment_char  :
						if not ( ( sep == None ) | ( sep == " " ) ) :
							line_elements = line.split( None ) #annotations are split with spaces
						if len( line_elements ) > 1 :
							last_character = len( line_elements[ 1 ] ) - 1
							if line_elements[ 1 ][ last_character ] == ":" :
								annotation_name = line_elements[ 1 ][ 0 : last_character  ]
								annotation = str("")
								for i in range( 2 , len( line_elements ) ):
									if i != len( line_elements ) - 1 :
										annotation = annotation + str(line_elements[ i ]) + " "
									else :
										annotation = annotation + str(line_elements[ i ]) #the last element has not space following
								self.annotations( annotation_name , annotation ) 
		if 'frames' in output.keys():
			try:
				self.input_values('frames',output['frames'])
			except:
				raise AttributeError('.load_data: chronological disorder in trajectory "' + file_name + '".')
		if 't' in output.keys():
			if 't_unit' not in self._annotations.keys():
				self.input_values('t',output['t'])
			else :
				self.input_values('t',output['t'],unit=self._annotations['t_unit'])
			self.input_values('t',output['t'])
		for item in output.keys():
			if (item != 'frames') | (item != 't'):
				if item == 'coord' :
					if 'coord_unit' not in self._annotations.keys():
						self.input_values(item,output[item])
					else :
						self.input_values(item,output[item],unit=self._annotations['coord_unit'])
				else :
					self.input_values(item,output[item])
		for a in [a for a in attrs.keys() if '_'+a not in self.__slots__]:
			if a not in self._annotations.keys():
				self.annotations(a,attrs[a])	
			else :
				raise AttributeError(a+' has been already annotated as: '+self._annotations[a])

	def fill(self):
		"""
		fill() fills attributes of missing frames with Nan
		"""
		non_empty_attributes = self.attributes()
		if 'frames' in non_empty_attributes: #Are frames empty?
			#Check if there are missing frames
			frame_intervals = self._frames[1:]-self._frames[0:(len(self._frames)-1)]
			if max(frame_intervals) > 1:
				for i in range(len(frame_intervals)-1,-1,-1):
					while frame_intervals[i] > 1:
						for attribute in non_empty_attributes:
							if attribute == 'frames' : #insert the missing frame
								x = insert(getattr(self,'_frames'),i+1,self._frames[i+1]-1)
								setattr(self,"_frames",array(x,dtype='int64')) # add the coords
							#insert a NaN in the attribute
							elif attribute == 't' : #insert the missing time point
								delta_t = min(self._t[1:]-self._t[0:(len(self._t)-1)])
								x = insert(getattr(self,'_t'),i+1,self._t[i+1]-delta_t)
								setattr(self,"_t",array(x,dtype='float64')) # add the coords
							elif attribute == 'coord' :
								x = getattr(self,'_coord')
								setattr(self,"_coord",array([\
										insert(x[0,],i+1,NaN),
										insert(x[1,],i+1,NaN)\
												],dtype='float64'))
							elif attribute == 'coord_err' :
								x = getattr(self,'_coord_err')
								setattr(self,'_coord_err',array([\
										insert(x[0,],i+1,NaN),
										insert(x[1,],i+1,NaN)\
												],dtype='float64'))
							else: 
								x = insert(getattr(self,'_'+attribute),i+1,NaN)
								setattr(self,"_"+attribute,array(x,dtype='float64')) # add the coords
						frame_intervals[i] = frame_intervals[i]-1 #The number of missing frames has been reduced by one
		elif 't' in non_empty_attributes: #Are times empty?
			#Check if there are missing frames
			time_intervals = self._t[1:]-self._t[0:(len(self._t)-1)]
			delta_t = float( self.annotations()[ 'delta_t' ] ) #old version:# min(time_intervals)
			time_intervals = time_intervals/delta_t

			if max(time_intervals) >= 2: #before (211101) it was > 1 insteaad of >= 2
				for i in range(len(time_intervals)-1,-1,-1):
					while time_intervals[i] >= 2:
						for attribute in non_empty_attributes:
							if attribute == 't' : #insert the missing time point
								x = insert(getattr(self,'_t'),i+1,self._t[i+1]-delta_t)
								setattr(self,"_t",array(x,dtype='float64')) # add the coords
							#insert a NaN in the attribute
							elif attribute == 'coord' :
								x = getattr(self,'_coord')
								setattr(self,"_coord",array([\
										insert(x[0,],i+1,NaN),
										insert(x[1,],i+1,NaN)\
												],dtype='float64'))
							elif attribute == 'coord_err' :
								x = getattr(self,'_coord_err')
								setattr(self,'_coord_err',array([\
										insert(x[0,],i+1,NaN),
										insert(x[1,],i+1,NaN)\
												],dtype='float64'))
							else: 
								x = insert(getattr(self,'_'+attribute),i+1,NaN)
								setattr(self,"_"+attribute,array(x,dtype='float64')) # add the coords
						time_intervals[i] = time_intervals[i]-1 #The number of missing frames has been reduced by one
		else:
			pass

	def attributes(self):
		"""
		attributes() reports the attributes of the trajectory that are
		not empty
		"""
		
		non_empty_attributes = []
		for s in self.__slots__[1:]:
			x = getattr(self,s)
			if (x.shape[x.ndim-1] > 0):
				non_empty_attributes.append(s[1:])
		return non_empty_attributes

	def annotations(self,annotation=None,string=''):
		if (annotation == None ) & ( not string ) :
			return self.__dict__()
		elif (annotation == None ) & ( not ( not string ) ) :
			raise AttributeError('You annotate something to no dictionary key!')
		elif (annotation != None ) & ( not string ) :
			if type( annotation ) == type( dict() ) :
				for key in annotation.keys() :
					self._annotations[ key ] = annotation[ key ]
			else :
				self._annotations[annotation] = string
		else:
			self._annotations[annotation] = string
	
	def assign_datasetID( self , path , pattern = 'Traj' ) :

		datasets =  [ f for f in os.listdir( path ) if pattern in f ] #list all the files in path that have pattern

		found_dataset = False

		for d in datasets :

			with open( path + '/' + d , 'r' ) as f :

				for line in f :

					t = re.search( str( self.frames()[ 0 ] ) + '.+' + str( self.coord()[ 0 ][ 0 ] ) + '.+' + str( self.coord()[ 1 ][ 0 ] ) , line )
					
					if t :
						
						found_dataset = True
						self.annotations()[ 'dataset' ] = d
						break

			f.close()

		# if the loop did not break by now, there was a problem 
		if not found_dataset :
			
			raise TypeError( 'assign_datasetID has not found the trajectory in any dataset. Check that the  path is correct, or that the frame() and coord() defined in the trajectory match the definitions in the dataset.\n assign_dataseID was searching for the string: ' + str( self.frames()[ 0 ] ) + '.+' + str( self.coord()[ 0 ][ 0 ] ) + '.+' + str( self.coord()[ 1 ][ 0 ] ) )

	def integral( self , what , scale = 1 , two_dimentional = False ) :

		x = getattr( self , '_'+what )
		
		xx = []

		if x.ndim == 1 :
			x0 = 0 # reset the intergal starting value to 0
			for i in range( 0 , len( x ) ) :
				if x[ i ] == x[ i ] : 
					xx.append( x0 + x[ i ] * scale )
					x0 = xx[ -1 ]
				else : 
					xx.append( NaN )
		else :
			if two_dimentional :
				x0 = 0 # reset the intergal starting value to 0
				for i in range( 0 , len( x[ 0 ] ) ) :

					if ( ( x[ 0 ][ i ] == x[ 0 ][ i ] ) & ( x[ 1 ][ i ] == x[ 1 ][ i ] ) ) : 
						xx.append( x0 + sign( x[ 0 ][ i ] ) * sqrt( x[ 0 ][ i ] ** 2 + x[ 1 ][ i ] ** 2 ) * scale )
						x0 = xx[ -1 ]
				
					else : 

						xx.append( NaN )
			else :
				for j in range( 0 , x.ndim ) :
					xx.append( [] )
					x0 = 0 # reset the intergal starting value to 0
					for i in range( 0 , len( x[ j ] ) ) :
						if x[ j ][ i ] == x[ j ][ i ] : 
							xx[ j ].append( x0 + x[ j ][ i ] * scale )
							x0 = xx[ j ][ -1 ]
						else : 
							xx[ j ].append( NaN )
		return( xx )
