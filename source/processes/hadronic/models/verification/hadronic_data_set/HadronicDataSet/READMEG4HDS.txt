                                The Geant4 Hadronic Data Set: G4HDS

Simone.Gilardoni@cern.ch
Vladimir.Grichine@cern.ch

Last Update: 20/11/03

1	Introduction

The Hadronic Data Set for Geant4 (G4HDS in the following) is a structure implemented 
to store and to manage data such as cross sections of hadronic interaction


2	Data Format

A data format for the G4HDS has been chosen as a compromise between different 
existing databases (like EXFOR or ENDF/B-VI).


2.1	G4HDS directory structure

The G4HDS format of the directory structure follows the rule:

/G4HADATASET/primary particle name/physics process/secondary particle name

where:
o	G4HADATASET is the environmental variable to be set before the first use of the 
G4HDS, which contains the G4HDS root path. The set of this variable is 
mandatory. Usually is $G4INSTALL/HadronicDataSet.
o	primary particle name is the name of the incident particle and it has to be 
declared as a G4string. The name coding follows the Geant4 particle name table. 
Example: G4string proton("proton");
o	physics process is the physics process whose data are stored inside the 
subdirectory structure. The list of the current permitted process can be found in 
the file $G4HADATASET/G4HDSProcesses.txt and at the end of this file. 
o	secondary particle is the name of the secondary particle produced by the physics 
process. It has to be declared as the primary particle name. In case that the 
process has no secondaries, like for the total reaction cross section, the secondary 
particle has to be declared as an empty string. 
Example: G4String nothing(""); 

The directory structure is copied from the G4HDS installation. 
It is possible to add to the structure any primary particle or physics process or secondary 
particle creating the proper directory inside G4HDS. This can be done if the 
environmental variable G4HADWRITTINGRIGHTS is set and equals 1.

Examples: 
/HadronicDataSet/proton/doublediff/neutron
/HadronicDataSet/proton/doublediff/proton 


2.2	File name format

Every file of G4HDS contains different process information defined by its path for a 
specific target defined by its filename. 
The target material has to be built as a G4Material or a G4Element or a G4Isotope.
This information is then processed to decode the filename.
The filename format is the following:
PhysicsprocessidZxxxxAxxxTargetid
where:
o	Physicsprocessid is the abbreviation for the physics process (dd for example stays 
for double differential cross section). See the file 
$G4HADATASET/G4HDSProcesses.txt for the list of known processes and their 
abbreviation.
o	Zxxxx is the effective atomic number Z of the G4Material or the G4Element or the 
Z of the G4Isotope.
o	Axxx is the effective nuclear weight A of the G4Material or the G4Element or the 
A of the G4Isotope.
o	Targeteid specifies if the file is associate to a G4Material (Targetid=mat) or to an 
G4Element (Targetid=ele) or to a G4Isotope(Targetid=iso)
 
Example:
Aluminium double differential cross section for impinging protons producing neutrons. 
 
density = 2.700*g/cm3; 
a = 26.98*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13, a, density);
  
filename = ddZ13.000A26.980mat

Full path:
/HadronicDataSet/proton/doublediff/neutron/ddZ13.000A26.980mat



2.3	Data file format

The format of the data contained into the file has been chose as follows according to the 
different physics process. 
Every file contains a header to define the source of the data (like articles, databases or 
collaborations). The header lines start with a #, which are considered as comments.
Special words, listed in the table below, identify different information blocks.
The data starts with the primary kinetic energy, identified by the special word Tkin, 
then the no word defines the number of bin for the information contained in the data.
The dimensions of the data are identified by special words.

Special Word		   Information		          Example
//                      Start comment line		  // this is a comment
#			Start comment line		  # this is a comment
Tkin			Primary kinetic energy	          Tkin	 1.13e+2
no			Number of bins			  no 4
angleUnit		Unit for angles			  angleUnit degree
energyUnit		Unit for energy			  energyUnit MeV
angleNo			Bin angle considered		  AngleNo 80
ddXscUnit		Unit for double diff. XS	  ddXscUnit barn/sr/MeV


A special comment line identifies the end of a data block inside the file:
#End data file

The new data relative to the same material and the same process but for different primary 
particle energy are appended after this line.



3	Functionalities

3.1	Writing function


G4HDS has been conceived to leave to the user the freedom to introduce data in any 
format and then translate them into the G4HDS format described above.
The virtual class G4VHadDataWritting has to contain the methods to write the data. 

3.2	Retrieving function

The retrieving methods are contained in the classes G4ReadHadprocess where process 
stays for one of the following:

o	DiffAngleXSC
o	DiffMomentumCXSC
o	DoubleDiffXSC
o	DiffMomentumCXSC
o	Multiplicity
o	InelasticXSC
o	TotalXSC
o	ElasticXSC

The scope of the methods is the identification of the data file. Then the data are retrieved 
by the method of G4HadDataReading class, which fills G4PhysicsTable with the 
recovered information.  

3.3	Searching function

The file searching function is implemented inside the class G4HadFileFinder.
The method to query for a file needs as parameter the starting searching path and the 
filename with the searching conditions.  
The path argument is used to fix the starting point of the search and to select either the 
primary particle or the process or the secondary particle. Setting path = G4HADATASET 
request a query in the entire G4HDS.
The filename argument is used to search for a specific file whose name has to be 
conforming to the rules explained in the data format section.  
The filename argument can contain selection characters to control the search, according 
to the following rules:
o	A * in any position means no selection on the specific field. 
Accepted filename: *Z13.00A26.980mat, ddZ*A26.980mat, ddZ13.00A*mat

o	A < or > after the Z character or the A character will select files according to the 
smaller than or greater than rule. 
Accepter filename: ddZ>13A26.980mat, ddZ<13A26.980mat, 
ddZ13.00A>26.980mat, ddZ13.00A<26.980mat

The query results are the full filename path, if any file is found inside the data set, which 
is written on the screen. 

Example:
 
G4String  temp;
G4String  path =  getenv("G4HADATASET");

temp = "*Z13.000A26.980mat";
G4HadFileFinder* findit=new G4HadFileFinder(path,temp)



4	Processes

Abbreviation		Directory name				Process description

dd 			doublediff				double differential cross section 
df			diff					differential cross section
el			elastic					elastic cross section
in			inelastic				inelastic cross section
mu			multeplicity				multeplicity distributions
to			total 					total reaction cross section
