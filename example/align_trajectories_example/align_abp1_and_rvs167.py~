from trajalign.align import align
from trajalign.average import load_directory

rvs167_trajectories = load_directory(
		path = 'abp1_and_rvs167' , 
		pattern = '.rvs167_data.txt$' ,
		comment_char = '%' , 
		dt = 0.2715 , 
		t_unit = 's' , 
		coord_unit = 'pxl' , 
		frames = 0 , 
		coord = ( 1 , 2 ) , 
		f = 3 , 
		protein = 'Rvs167-GFP' , 
		date = '10/10/13' , 
		notes = 'the trajectory of the target protein')

abp1_trajectories = load_directory(
		path = 'abp1_and_rvs167' , 
		pattern = '.abp1_data.txt$' ,
		comment_char = '%' , 
		dt = 0.2715 , 
		t_unit = 's' , 
		coord_unit = 'pxl' , 
		frames = 0 , 
		coord = ( 1 , 2 ) , 
		f = 3 , 
		protein = 'Abp1-mCherry' , 
		date = '10/10/13' , 
		notes = 'the trajectory of the reference protein')

align( path_target = 'rvs167.txt' , path_reference = 'abp1.txt' , ch1 = rvs167_trajectories , ch2 = abp1_trajectories )
