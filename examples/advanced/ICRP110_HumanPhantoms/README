     =======================================================================
                    Geant4 - ICRP110_HumanPhantoms Example 
     =======================================================================

The ICRP110_HumanPhantoms example is developed and mantained by Susanna Guatelli, Matthew Large and Alessandra Malaroda,
Centre For Medical Radiation Physics (CMRP), University of Wollongong, NSW, Australia, and John Allison, Geant4 Associates International 
and University of Manchester, UK.

  Contacts:     
    - susanna@uow.edu.au 
    - mjl970@uowmail.edu.au
    - malaroda@uow.edu.au
    - John.Allison@g4ai.org

The example is based on the extended/medical/DICOM example

The authors acknowledge that this application of the ICRP110 human phantoms have been implemented in Geant4 with the kind permission of 
the International Commission on Radiological Protection (ICRP). 

----------------------------------------------------------------------------------------------------
--------------------------------------> Introduction <----------------------------------------------
----------------------------------------------------------------------------------------------------

This application models the ICRP110 reference computational human phantoms [1] in a Geant4 simulation and calculates
the dose in individual voxels and in entire organs. 

The human male phantom, provided kindly by the ICRP, is created from a whole-body clinical CT image set of a 38yr old 
individual with height 176 cm and mass approximately 70 kg. Similarly, the human female phantom was created from a set of 
whole body CT images of a 43yr old individual with height 163 cm and weight 60 kg. The CT scans were acquired with both 
individuals laying supine and with arms resting parallel alongside the body. Both sets of CT data were then scaled to 
closely approximate the ICRP adult Reference Male and Reference Female, defined in previous ICRP publications [2, 3]. 
 
	[1] HG Menzel, C Clement, and P DeLuca. ICRP publication 110. "Realistic reference phantoms:
	an icrp/icru joint effort: A report of adult reference computational phantoms", Annals of the
	ICRP, 39(2):1, 2009. URL: http://www.icrp.org/publication.asp?id=icrp%20publication%20110.
 
  [2] Valetin J 2002 Basic anatomical and physiological data for use in radiological protection: 
  reference values: ICRP Publication 89 Ann. ICRP vol. 32 (Oxford: Elsevier) pp 1-277.

  [3] Valetin J 2007 The 2007 recommendations of the international commission on radiological 
  protection Ann. ICRP vol 37 (Oxford: Elsevier) pp 1-133.

The table below summarises the key features of the male and female voxelised human phantoms.
    
	PROPERTY			       	AM		AF
	_____________________________________
	Height (m)			    	1.76	1.63	
	
	Mass(Kg)				      73.0	60.0	
	
	Slice Thickness(mm)		8.0		4.84
	
	Voxel in-plane-			  2.137	1.775
	   -resolution (mm)				
	   
	Voxels along x 			  254		299		
		(i.e. columns)
		  
	Voxels along y			  127		137
		(i.e. rows)
		
	Number of Slices		  222		348
		(i.e. along z)
	______________________________________	

----------------------------------------------------------------------------------------------------
------------------------------> Application Sub-Folder Structure <----------------------------------
----------------------------------------------------------------------------------------------------

 - '/src': where the source .cc files are stored

 - '/include': where header .hh files are stored

 - '/ICRPdata': where the phantom data files (*.dat) and slice files are stored.
 It is downloaded automatically from URL https://cern.ch/geant4-data/datasets/examples/advanced/ICRP110Phantoms/ICRPdata.tar.gz
 during the configuration via cmake.
 
 Phantom data files containing the voxelisation of each phantom, as well as files 
 containing the definitions of the phantom organs and materials used within geant4 
 code can be found in the folder /ICRPdata. 
 
 All data files used for this phantom were obtained from the ICRP's website on publication 110 under "Supplementary Data"  
        - https://www.icrp.org/publication.asp?id=ICRP%20Publication%20110. 
        
----------------------------------------------------------------------------------------------------
----------------------------------> ICRP110Phantoms Data <------------------------------------------
----------------------------------------------------------------------------------------------------

Within the '/ICRPdata' directory, the following sub-directories are contained:

        -> /ICRPdata/                  : contains '*Data.dat' files which list the number of phantom slices to 
                                         simulate and the order in which to stack the phantom slices.
                                         
        -> /ICRPdata/ICRP110_g4dat/AM/ : contains the individual male phantom slice files. 
                                         
        -> /ICRPdata/ICRP110_g4dat/AF/ : contains the individual female phantom slice files.
                                         
        -> /ICRPdata/ICRP110_g4dat/P110_data_V1.2 

The final directory contains the raw ICRP110 phantom data as obtained from the ICRP110 publication website [1]; 
5 files within folders for the AM and AF phantoms are given. These files are described as follows in the 
supplementary data's included README file. 

 The array of organ identification numbers (in ASCII format); the file names are: 
		AM.dat
		AF.dat

 A list of individually segmented structures, their identification numbers, and assigned media (Appendix A in ICRP110); the file names are: 
		AM_organs.dat
		AF_organs.dat
   
	A list of the media, their elemental compositions and densities (Appendix B in ICRP110); 
  the file names are: 
		AM_media.dat
		AF_media.dat
  
  The mass ratios of bone constituents (trabecular bone, red and yellow bone marrow) in the spongiosa regions; 
  the file names are: 
		AM_spongiosa.dat
		AF_spongiosa.dat
	
  The mass ratios of blood in various body tissues; the file names are: 
		AM_blood.dat
		AF_blood.dat

The primary data files AM.dat and AF.dat contain an array of organ identification numbers ranging from 0 to 141.
Each number respresents the organ associated with each voxel within the phantom. Within these files, the organ IDs 
are listed slice by slice, within each slice row by row, within each row column by column. That means, the column 
index changes fastest, then the row index, then the slice index - in other words, the phantom voxels first increase 
along x, then along y and finally along z. Slice numbers increase from the toes up to the vertex of the body; 
row numbers increase from front to back; and column numbers increase from right to left side.
 
For use in this application, the original AM.dat and AF.dat files containing the organ identification numbers of 
all voxels of the phantom were sub-divided into many files with each representing a single phantom slice along z. 
As such, each file represents a 2D phantom slice containing x,y voxel positions and organ identification numbers 
of each voxel. This allows for subsections of the phantom to be simulated as required by the user, removing the 
need to simulate the entire phantom every time when this may not nessecrily be needed by the user. This also will 
allow for reductions in the simulation time depending on what portion of the total phantom is simulated by the user. 
This feature was achieved via a code developed by Dr Alessandra Malaroda, University of Wollongong, Australia in 2017. 

The AM human phantom is voxelised in x,y,z with 254 x 127 x 222 voxels with dimensions 2.137 x 2.137 x 8 mm.
The AF human phantom is voxelised in x,y,z with 299 x 137 x 348 voxels with dimensions 1.775 x 1.775 x 4.84 mm.

----------------------------------------------------------------------------------------------------
---------------------------------------> How to compile and run <-----------------------------------
----------------------------------------------------------------------------------------------------

- Create a build folder for the phantom run
		  % mkdir build/
		  
- Navigate to inside the build folder and initialise Geant4 
		  % cmake ../

  The ICRP110 phantom data will be automatically downloaded from https://cern.ch/geant4-data/datasets/examples/advanced/ICRP110Phantoms/ICRPdata.tar.gz
  
- Compile and link to generate the executable (in your CMAKE build directory):
 	      % make
  This should make two executables - ICRP110phantoms and ICRP110standalone.

- Execute the application in 'interactive' mode with visualization:
 	      % ./ICRP110phantoms
              
- Execute the "standalone" application in 'interactive' mode with visualization:
 	      % ./ICRP110standalone
  This allows you to visualise the phantom without the overhead of the run manager and initialising all the physics tables.
  Of course, you cannot run or visualise trajectories.
              
- Execute the application in 'batch' mode from macro files:
        % ./ICRP110phantoms female_head.in

-----------------------------
    AVAILABLE MACRO FILES 																																	
-----------------------------	
For the users convenience, macro files have been created which are designed to construct partial head 
and trunk phantoms for both the male and female models. These macro files can be called upon in batch
mode when executing the application as specified above. If the user wishes to construct a completed/full
male or female phantom, the macros male.in and female.in can be called upon, respectively. 

 - male_head.in/female_head.in   : Creates a partial head phantom for the male and female, respectively.
 - male_trunk.in/female_trunk.in : Creates a partial trunk phantom for the male and female, respectively.  
 - male.in                       : Creates full male ICRP110 phantom. This can be modified along with 'ICRPdata/MaleData.dat'
                                   if the user wishes to create their own custom partial phantom section.
 - female.in                     : Creates full female ICRP110 phantom. This can be modified along with 
                                  'ICRPdata/FemaleData.dat' if the user wishes to create their own custom partial phantom section.  
 - openGLVis.mac                 : macro for visualisation with openGL. 
 - vis.mac (default)             : Executed by default when the simulation is run in 'interactive' mode.	
 - primary.mac                   : Contains the definition of the primary radiation field.

At the very top of the various '.in' macro files (pre-initialization), there are a series of commands 
which define the sex and section of the phantom to create. These commands are listed below:

  o /phantom/setPhantomSex <option> : Passes sex of phantom to Detector Construction
  o /phantom/setScoreWriterSex <option> : Passes sex of phantom to User Score Writer
  
  o /phantom/setPhantomSection <option> : Passes section of phantom to Detector Construction
  o /phantom/setScoreWriterSection <option>  Passes section of phantom to User Score Writer
  
Available options for the first 2 commands are: male or female.
Avalable options for the last 2 commands are: head, trunk or full.  
  
In the event that the macro called upon by the user when executing the application in 'batch' mode
does not contain these commands (default case), the application sets phantom sex to female and the section as the head. 

WARNING: the phantom model can be chosen only in the initialization phase of the simulation!!!
It cannot be changed during the run session. This feature will be implemented in the next future. 

----------------------------------------------------------------------------------------------------
----------------------------------> Creating a Custom Phantom <------------------------------------
----------------------------------------------------------------------------------------------------

If the user wishes to construct a customised section of the phantom (i.e. a single slice, the legs, etc),
he/she has to create a specific macro or edit the ones provided. The recommended method for a custom male 
phantom is outlined as follows. 

The user should edit the macro 'male.in' and the data file 
'MaleData.dat'. Firstly, in 'FemaleData.dat', there are 2 simple ways in which the user can
select a custom range of phantom slices to simulate: 

1. The very first entry of each Data.dat indicates how many slices to simulate. 
   Changing this number will determine the number of slices to construct.
   
2. Further down in the Data.dat files (beginning at line 61) is the name of the first slice to simulate, followed 
   by successive slices. Changing the slice file orders here will allow various subsections of the human 
   phantom to be simulated. As an indication the following phantom subsections have been identified for the  
   male phantom below.
   
		--> AM_Slice1.g4dat to AM_Slice20.g4dat: Feet to ankles
		
		--> AM_Slice21.g4dat to AM_Slice121.g4dat: Ankles to hips
		
		--> AM_Slice169.g4dat: Single chest slice with good visualisation 
                           of lungs, ribs, heart.
		
		--> AM_Slice182.g4dat to AM_Slice222.g4dat: Neck and Head
		
       NOTE: o Always order phantom slices beginning with the lowest number and increasing 
               in slice number going down the .dat files.
    	       o Always use consecutive/adjacent slices when simulating multiple slices.
             o The default number of slices for both male and female phantoms is set to 10
               and starts at the feet of each phantom.	
 
Once the user customises the MaleData.dat/FemaleData.dat (for example starting from the full phantoms macros), 
he/she has also to fix appropriately the scoring mesh in male.in/female.in.

----------------------------------------------------------------------------------------------------
------------------------------> Scoring Mesh and the User Score Writer <----------------------------
----------------------------------------------------------------------------------------------------

The macro primary.mac defines the radiation beam type, energy, direction and geometry. The UI commands of the
General Particle Source should be used to change the radiation field. The macros male.in and female.in contain 
the /run/beamOn command and can call upon the radiation beam definition through the UI command
'/control/execute primary.mac'. 

Within male.in and female.in, a scoring mesh is defined which records the dose deposition within each individual 
phantom voxel. The size of the scoring mesh is defined in line 54 of the male.in/female.in files, and must be defined 
to match the constructed phantom dimensions (whole or partial) defined in the according '/ICRPdata/*Data.dat' file. 

The mesh dimensions are defined as half-dimensions in x,y,z - meaning a defined scoring mesh x-dimension of 100mm will construct 
a scoring mesh spanning from -100mm to +100mm in the geometrical world in which the phantom lies. Furthermore, for the completed 
male phantom which has dimensions along x,y,z of 542.798 x 271.399 x 1776 mm, the scoring mesh half-dimensions should be defined 
as 271.399 x 135.6995 x 888. mm. The number of bins or divisions to segment the mesh into is then defined in line 51. These 
should match the number of phantom voxels in x,y,z which are defined in the MaleData.dat and FemaleData.dat files in the '/ICRPdata' 
directory.  

If the user edits the MaleData.dat or FemaleData.dat files to change the number of z-slices simulated in a run, they must also edit 
the scoring mesh dimensions and number of bins to ensure it correctly scores their defined phantom. To do so, the user will typically 
only have to edit lines 54 and 55 of the male.in or female.in macro files.

After completion of a simulation run, the phantom mesh records the deposited dose in each phantom voxel and outputs the data to a text file named 
"PhantomMesh_Dose.txt". This text file lists the x,y,z positional number of the voxel in the phantom and the dose recorded within that voxel (in Gy). 

The output PhantomMesh_Dose.txt file is created by the User Score Writer class defined in the source code ICRP110UserScoreWriter.cc. In the same class the dose 
in the voxels is analysed and associated to organs. 

A final output file "ICRP.out" is then created which contains the total dose delivered to each organ.   

----------------------------------------------------------------------------------------------------
----------------------------------------> Further Info <--------------------------------------------
----------------------------------------------------------------------------------------------------

-------> ColourMap.dat <--------

This file located in the build directory assigns G4colours to the 53 phantom materials.
The user may edit these as they wish for visualistion purposes. 

----------> Physics <-----------

The QGSP_BIC_HP Physics List is adopted. The user may want to change the
cut of production of secondary particles. 

-----> Primary particles <------

The G4 General Particle Source (gps) is used to generate primary radiation field.
Macro primary.mac contains the definition of the primary radiation field.   
