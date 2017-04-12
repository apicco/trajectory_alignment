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
from numpy import cos
from numpy import insert
from numpy import NaN
from numpy import insert
from numpy import nanmean
from numpy import nanmax
from numpy import nanmin
from numpy import float64 
from numpy import convolve
from math import isclose
import copy as cp

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

		.t(), .coord(), .f() and .n() have related attributes containing the 
		measured errors, if known. The errors are called with .t_err(), 
		.coord_err(), .f_err() and .mol_err() respectively
		
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
		"f", "n", "mol","t_err", "x_err", "f_err", "mol_err".

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
	
	__slots__ = ['_annotations','_frames','_t','_coord','_f','_mol','_n','_t_err','_coord_err','_f_err','_mol_err']
	 

	def __init__(self,**annotations):

		#Trajectory main attributes 
		self._annotations = annotations
		self._frames = array([],dtype='int64')
		self._t = array([],dtype='float64')
		self._coord = array([array([]),array([])],dtype='float64') #equivalent to matrix in python3.5 provided I use @ operator
		self._f = array([],dtype='float64')
		self._mol = array([],dtype='float64')
		self._n = array([],dtype='float64')

		#Trajectory error attributes
		self._t_err = array([],dtype='float64')
		self._coord_err = array([array([]),array([])],dtype='float64') #equivalent to matrix in python3.5 provided I use @ operator
		self._f_err = array([],dtype='float64')
		self._mol_err = array([],dtype='float64')


	def __dict__(self):
		return self._annotations

	def __len__(self): #the number of timepoints in the trajectory
		if (len(self._t) > 0): return len(self._t)
		else: return len(self._frames)

	def __repr__(self):
		if (len(self)) == 0 :
			output = 'The trajectory is empty!\n'
			if (len(self._annotations)):
				output += '#' + '-' * len(output) + '\n'
				for name, item in self._annotations.items():
					output += '# ' + name + ': ' + item + '\n'
			return output
		else :
			table = []
			names = []
			output = '#'
			for s in self.__slots__[1:]:
				x = getattr(self,s)
				if (x.shape[x.ndim-1] > 0):
					if s in ('_coord','_coord_err'):
						#x coord
						table.append(x[0])
						if len(self._annotations['coord_unit']) > 0:
							names.append('x' + s[6:] + ' (' + self._annotations['coord_unit'] + ')')
						else:
							names.append('x' + s[6:])
						#y coord	
						table.append(x[1])
						if len(self._annotations['coord_unit']) > 0:
							names.append('y' + s[6:] + ' (' + self._annotations['coord_unit'] + ')')
						else:
							names.append('y' + s[6:])
					else: 
						table.append(x)
						if s[1:] + '_unit' in self._annotations.keys():
							if len(self._annotations[ s[1:] + '_unit' ]) > 0:
								names.append( s[1:] + ' (' + self._annotations[s[1:] + '_unit' ] + ')' )
							else: 
								names.append(s[1:])
						else:
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
			for row in transpose(table):
				output_row = ''
				for element in row:
					output_row += str(element).rjust(col_width)
				output += output_row+'\n'
			#print annotations
			output += '#'+'-'*(len(names)*col_width-1)+'\n'#nice separator
			for name, item in self._annotations.items():
				output += '# '+name+': '+item+'\n'
			return output


	#Getters

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
				if a is not 'range':
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
					if a == 't_err' : 
						output.input_values(a,self.t_err(new_items))
					if a == 'coord_err' : 
						output.input_values(a,self.coord_err(new_items))
					if a == 'f_err' : 
						output.input_values(a,self.f_err(new_items))
					if a == 'mol_err' : 
						output.input_values(a,self.mol_err(new_items)) 
			except IndexError:
				print('Indexes in range are out of bounds')

		return output

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
			elif ((name=='t') & (len(self._frames)==0)):
				setattr(self,"_"+name,array(x,dtype='float64')) # add the time; no frames present yet
				if ('t_unit' in self._annotations.keys()):
					if (self._annotations['t_unit'] == ''): self._annotations['t_unit'] = unit
				else: self._annotations['t_unit'] = unit
			elif (name=='coord'):
				if (len(x[0])==(len(self._frames) | len(self._t))):
					setattr(self,"_"+name,array(x,dtype='float64')) # add the coords
					if ('coord_unit' in self._annotations.keys()):
						if (self._annotations['coord_unit'] == ''): self._annotations['coord_unit'] = unit
					else: self._annotations['coord_unit'] = unit
				else:
					raise AttributeError('The input length of the attribute does not match frames or time array of non-zero length')
			elif (name=='coord_err'):
				if (len(x[0])==(len(self._frames) | len(self._t))):
					setattr(self,"_"+name,array(x,dtype='float64')) # add the coords
					if ('coord_unit' in self._annotations.keys()):
						if (self._annotations['coord_unit'] == ''): self._annotations['coord_unit'] = unit
					else: self._annotations['coord_unit'] = unit
				else:
					raise AttributeError('The input length of the attribute does not match frames or time array of non-zero length')
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
		self._f = ( self._f - nanmin(self._f) ) / ( nanmax(self._f) - nanmin(self._f) )
		#if the attribute _f_err is not empty, then propagate the errors accordingly 
		if ( self._f_err.shape[ self._f_err.ndim - 1 ] > 0 ) :
			self._f_err = self._f_err / ( nanmax(self._f) - nanmin(self._f) )
			#the aim of this function is only to rescale the fluorescence intensity
			#between 0 and 1, hence we do no propagate the error of the max(self._f)
			#and min(self._f).

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
		center_mass(): centers the trajectory on its center of mass
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
					( ( t  > self._t[0] ) | isclose( t , self._t[0] ) ) & 
					( ( t < self._t[len(self)-1] ) |  isclose( t , self._t[0] ) )
						):
				new_t = array([ i for i in self._t if ( (i > t) | isclose(i,t) ) ])
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

				number_of_new_frames = float64( (self._t[0]  - t)/delta_t )
				
				#the difference between the time at which the trajectory starts and the time
				#of the new starting point is divided by the time interval delta_t to compute the
				#number of correspondent frames, which is a float. Because of the representation of 
				#float numbers a time interval which corresponds, for exampe to an exact number of 
				# 5 frames might be computed as  4.999999999999999 instead of 5, which would round 
				#to the wrong integer. Here, we makes sure that the number_of_new_frames is 
				#rounded to the right integer.
				if int( number_of_new_frames + 0.01 ) != int( number_of_new_frames) :
					number_of_new_frames = int( number_of_new_frames + 0.01)
				else :
					number_of_new_frames = int( number_of_new_frames)
				new_t = [ self._t[ 0 ] - i * delta_t for i in range( number_of_new_frames , 0 , -1 ) ] 
				self._t = insert(self._t,0,new_t)
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
					( ( t  > self._t[0] ) | isclose( t , self._t[0] ) ) & 
					( ( t < self._t[len(self)-1] ) |  isclose( t , self._t[0] ) )
						):
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
				def float_range(x,y,step): #an inline function to compute the float range
					float_range_output = [];
					while ( (x < y) | isclose(x,y) ):
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
								raise TypeError('The comment_char might be ill-defined. Default is "#".')
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
			delta_t = min(time_intervals)
			time_intervals = time_intervals/delta_t

			if max(time_intervals) > 1:
				for i in range(len(time_intervals)-1,-1,-1):
					while time_intervals[i] > 1:
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
		if (annotation == None) & (string == '') :
			return self.__dict__()
		elif (annotation ==None) & (string != '') :
			raise AttributeError('You annotate a string to nothing!')
		elif (annotation !=None) & (string == '') :
			if type( annotation ) == type( dict() ) :
				for key in annotation.keys() :
					self._annotations[ key ] = annotation[ key ]
			else :
				return( self._annotations[annotation] )
		else:
			self._annotations[annotation] = string

