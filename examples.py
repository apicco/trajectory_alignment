#####################################
# examples of usage of the software:
# the trajectory class and its 
# functionalities
#####################################

from tracking.format.traj import Traj
#from tracking.align.align import align
from matplotlib import pyplot as plt

#read the content of the class Traj
dir(Traj)

#many functions! all those starting with __*__ are system variables of methods. #The use of the others can be accessed with:
print(Traj.__doc__)

#an empty trajectory is defined as
t = Traj()

#indeed, it is empty:
print(t)

#when a trajectory is defined, we can enter annotations in the form of a dictionary:
t = Traj(what="my first trajectory",mood="today is a beautiful day")
print(t)

#these annotations can be read
print(t.annotations()['what'])

#or modified
t.annotations()['mood']= 'all this seems pretty complicated. Depressing!'
print(t)

#we can fill some attributes of the trajectory, like the time:
t.input_values('t',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
print(t)

#note that there is a new attribute that has been added but is empty:
#the time unit (t_unit) . It is a bad habit to keep it empty, and 
#should be filled:
t.annotations()['t_unit']= 's'
print(t)

#let's add the spatial coordinates, and do not forget the units:
t.input_values('coord',[[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]],unit='apples')
print(t)

#let's load a real trajectory. We define the new trajectory first
d = Traj(experiment = 'Sla1-GFP', date = '25/12/2000', temperature = 'hot!')
d.load("example_trajectories/01.data",frames = 0, coord = (1,2), f = 3)

#let's correct the comment character, but I do not remember the variable name:
print(d.load.__doc__)
#comment_char!
d.load(file_name='example_trajectories/03.data',frames = 0, coord = (1,2), f = 3, comment_char='%', weather = 'clear and sunny')
print(d)

#we add the units
d.annotations()['coord_unit']= 'pxl'
print(d)

#let's compute the time from the frame numbers, we know that each frame is separated by 0.1045 s
d.time(0.1045,'s')

#there are some missing frames, easy! let's fill them:
d.fill()
print(d)

#let's save our trajectory, as a .txt file:
d.save('my_first_trajectory')

#let's rotate the trajectory (in radiants)
from numpy import pi
d.rotate(pi/2)

#let's shift it
d.translate((1000,-1000))

d.save('my_first_trajectory_rototraslated')
print(d)

#what are the first 3 time values:
print(d.t(0,1,2))
print(d.t([0,1,2]))
print(d.t(range(0,3)))

#what are the first 3 coordinate values:
print(d.coord(0,1,2))
print(d.coord([0,1,2]))
print(d.coord(range(0,3)))

#when the trajectory starts
print(d.start())

#when the trajectory ends
print(d.end())

#what is its length?
print(len(d))

#what is the trajectory lifetime?
print(d.lifetime())

#select only a part of the trajectory comprised in the range inputed
dshort = d.extract(10,11)
print(dshort)

dshort2 = d.extract(range(10,30))
#start dshort2 from t > 0
dshort2.start(0)
#end dshort2 at t < 36
dshort2.end(36)
print(dshort2)
