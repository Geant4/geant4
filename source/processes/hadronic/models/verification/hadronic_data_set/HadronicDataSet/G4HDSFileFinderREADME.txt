Geant4 Hadronic Set File Finder
Authors: S. Gilardoni, V. M.Grichine

1. Purpose of File Finder 

The Geant4 Hadronic Data Set (G4HDS) is a collection of files, which contain hadronic 
cross sections. The scope of the File Finder class is to search for a specified file inside the 
data set according to a specified query.

2. How to use File Finder

The file searching function is implemented inside the class G4HadFileFinder.
The method to query for a file needs as parameter the starting searching path and the 
filename with the searching conditions.  
The path argument is used to fix the starting point of the search and to select either the 
primary particle or the process or the secondary particle. Setting path = G4HADATASET 
request a query in the entire G4HDS.
The filename argument is used to search for a specific file whose name has to be 
conforming to the rules explained in the data format section of the user manual.  
The filename argument can contain selection characters to control the search, according 
to the following rules:
o	a * lists all the files contained in the G4HDS;
o	a * in any position means no selection on the specific field. 
Accepted filename: *Z13.000A26.980mat, ddZ*A26.980mat, 
ddZ13.000A*mat;
o	a < or > after the Z character or the A character will select files according to 
the smaller than or greater than rule; 
Accepter filename: ddZ>13A26.980mat, ddZ<13A26.980mat, 
ddZ13.000A>26.980mat, ddZ13.000A<26.980mat
o	the only mixed selections allowed are: 
no selection on the process and the Z or the A 
Accepted filename: *Z*A26.980mat, *Z13.000A*mat, ddZ*A*mat;

The query results are the full filename path, if any file is found inside the data set, which 
is written on the screen. 

Example:
 
char* temp;
char* path =  getenv("G4HADATASET");

temp = "*Z13.000A26.980mat";
G4HadFileFinder* findit=new G4HadFileFinder(path,temp);  

3. Limitation of the File Finder
o	The digits after value for the Z and the A have to be always three, even if they are 
zeros.
o	There is no search possibilities the target nature (material,nucleus, etc)

