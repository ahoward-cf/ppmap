# ppmap
Code for column density mapping with extra dimensions (temperature and dust opacity index).


PPMAP is still in heavy development. If you would like to run PPMAP, please contact K Marsh directly. 

Current Implimentation is for the Super Computing Wales Hawk architecture. 

===================================================================================================================================
HAWK Installation 

Clone repository to /home/[USERNAME]/
You should then have the following paths:

/home/[USERNAME]/ppmap/
/home/[USERNAME]/ppmap/source/

In /home/[USERNAME]/ppmap/ run 'sh make'. This should compile and place the ppmap and premap applications in /home/[USERNAME]/ppmap/.

Open /home/[USERNAME]/ppmap/run_ppmap.sh in your text editor. 
Replace instances of [USERNAME] with your HAWK login. 
Adjust the time required (header line #SBATCH -t d-hh:mm).

Open /home/[USERNAME]/ppmap/[FIELD]_premap.inp. Edit as needed. Rename to your field name. 

Open /home/[USERNAME]/ppmap/run_[FIELD]_all. Replace [USERNAME] with HAWK login. Replace [FIELD] with your field name. 
