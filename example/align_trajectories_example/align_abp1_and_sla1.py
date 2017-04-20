from trajalign.align import align
from trajalign.average import load_directory

sla1_trajectories = load_directory(
		path = 'abp1_and_sla1' , 
		pattern = '.sla1_data.txt$' ,
		comment_char = '%' , 
		dt = 0.2657 , 
		t_unit = 's' , 
		coord_unit = 'pxl' , 
		frames = 0 , 
		coord = ( 1 , 2 ) , 
		f = 3 , 
		protein = 'Sla1-GFP' , 
		date = '01/01/00' , 
		notes = 'the trajectory of the target protein')

abp1_trajectories = load_directory(
		path = 'abp1_and_sla1' , 
		pattern = '.abp1_data.txt$' ,
		comment_char = '%' , 
		dt = 0.2657 , 
		t_unit = 's' , 
		coord_unit = 'pxl' , 
		frames = 0 , 
		coord = ( 1 , 2 ) , 
		f = 3 , 
		protein = 'Abp1-mCherry' , 
		date = '01/01/00' , 
		notes = 'the trajectory of the reference protein')

align( path_target = 'sla1.txt' , path_reference = 'abp1.txt' , ch1 = sla1_trajectories , ch2 = abp1_trajectories , fimax2 = True )
