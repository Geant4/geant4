//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//


#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/dir.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>


#include <fstream>
#include <iomanip>
#include <iostream>
#include <strstream> 
#include <sstream>
#include "globals.hh"
#include "G4ios.hh"


#include "G4HadFileFinder.hh"

////////////////////////////////////////////////////////////////
//
//

G4HadFileFinder::G4HadFileFinder(G4String& currdir, G4String& G4filetofind)
{

  char filetofindone;
  char* pfiletofindone=&filetofindone;

  char directory;
  char* pdirectory=&directory;

  char nameChar[150] = {""};
 
  std::ostrstream ost(nameChar, 150, std::ios::out);
  ost << G4filetofind.c_str() ;

  pfiletofindone=nameChar;

  char nameDir[150] = {""};

  std::ostrstream ostdir(nameDir, 150, std::ios::out);
  ostdir << currdir.c_str();

  pdirectory=nameDir;
  
  char* workingresults = getenv("G4HADWORKING");
  G4int errnum = chdir(workingresults);
  
  
  std::ofstream filecheck("G4HadQuery.txt", std::ios::out);
  filecheck.close();
  
  // G4cout<<currdir<<nameChar<<G4endl;


  char* workingpath = getenv("G4HADATASET");
  if ( !workingpath )
  {
    G4String excep = "";
    excep += "G4HADATASET environment variable not set";
    G4Exception(excep);
  }
  errnum = chdir(pdirectory);

  G4cout << "File to find " << G4filetofind << " " << pdirectory << G4endl;

  Crawldir(pdirectory,G4filetofind);

  errnum = chdir(workingresults);
}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4HadFileFinder::~G4HadFileFinder()
{
}

////////////////////////////////////////////////////////////////
//
//

void G4HadFileFinder::Crawldir(char* currdir, G4String& filetofind)
{
  /* Open the directory */
  DIR *dir;
  struct dirent *de;
  struct stat st;
  char* totaldir;
  char tempname[MAXNAMLEN];


  G4String helpme="HELPME";

  // char* helpme="HELPME";

  
  if ( strcmp(filetofind,helpme) != 0 ) 
  {
    chdir (currdir);
    totaldir = getcwd(tempname,MAXNAMLEN);
    G4cout << totaldir << G4endl;  // what we see on screen each time, just for check
    dir = opendir(".");

    while ( (de = readdir(dir)) )
    {
      lstat( de->d_name, &st);

      if (S_ISDIR(st.st_mode)) 
      {
	if ((strcmp(de->d_name, ".") != 0) && (strcmp(de->d_name, "..") != 0)) 
        {
	      Crawldir(de->d_name, filetofind);
	}
      }
      else 
      {
	 G4String tempfilename(de->d_name);
	 G4StripFile(tempfilename, totaldir, filetofind);
      }
    }
    closedir(dir);
	//  G4cout<<'Switching back one dir level'<<G4endl;
	chdir("..");
  }
  else
  {
    char* workingpath = getenv("G4HADATASET");

    if ( !workingpath )
    {
       G4String excep = "G4HADATASET ";
       excep += "environment variable not set";
       G4Exception(excep);
    }
    //  G4int errnum = chdir(workingpath);

    std::ifstream fin("G4HDSFileFinderREADME.txt");
    std::filebuf* finbuf = fin.rdbuf();
	
    if ( !( finbuf->is_open() ) )
    {
      G4String excep = "G4HadDataSet -  G4HDSFileFinderREADME.txt not found. Did you delete something?";
      excep +="The help file G4HDSFileFinderREADME.txt is no more where it should be";
      G4Exception(excep);
    }
    else
    {
      G4String tmpString;

      while (!fin.eof())
      {
	getline(fin,tmpString);
	G4cout << tmpString << G4endl;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void G4HadFileFinder::G4StripFile(G4String& filename, char* dirposition, G4String& findmefile)
{
  G4String control = "*", pathfilename = "File not found", originalfile=filename;
  G4String pathfile(dirposition);
  G4String temp, process, A, Acondition, Z, Zcondition, target;
  G4String oprocess, oA, oAcondition, oZ, oZcondition, otarget;

  G4int zposition = -1, aposition = -1, targetposition = -1, replacelength = -1;  
  G4int processok = 0, Zok = 0, Aok = 0, targetok = 0, sumcheck = 0;

  G4double Zvalue = 0, Avalue = 0, oZvalue = 0, oAvalue = 0;


  char* workingresults = getenv("G4HADWORKING");
  G4int           errn = chdir(workingresults);
  G4String filenameout = "G4HadQuery.txt";

  // G4cout << filenameout << G4endl;

  std::ofstream fileout(filenameout, std::ios::app);  // app

  //  fileout.open(filenameout, std::ios::app);

  if ( (filename == findmefile) || (findmefile == control) )
  {
      G4cout << "File found in" << G4endl;
      G4cout<< dirposition << "/"<< filename << G4endl;
      fileout<< dirposition << "/"<< filename << std::endl;
  }
  else
  {

       // findmefile analysis

    zposition=findmefile.first('Z');
    process=findmefile.substr(0,zposition);
       
       //       G4cout << process << G4endl;

    zposition=findmefile.find('Z');
    aposition=findmefile.find('A');
       
       //G4cout << zposition << " " << aposition << G4endl;
       
    replacelength=(aposition-zposition)-1;
    zposition+=1;

    Z=findmefile.substr(zposition,replacelength);
       //       G4cout << Z << G4endl;
       
    targetposition=findmefile.length()-3;
    replacelength=(targetposition-aposition)-1;
    aposition+=1;

    A=findmefile.substr(aposition,replacelength);
       //       G4cout << A << G4endl;

    target=findmefile.substr(targetposition,3);
       //       G4cout << target << G4endl;


       // find from dir analysis

    zposition=filename.first('Z');
    oprocess=filename.substr(0,zposition);

       //       G4cout << oprocess << G4endl;

    zposition=filename.find('Z');
    aposition=filename.find('A');
       
    replacelength=(aposition-zposition)-1;
    zposition+=1;

    oZ=filename.substr(zposition,replacelength);
       //       G4cout << oZ << G4endl;
       
    targetposition=filename.length()-3;
    replacelength=(targetposition-aposition)-1;
    aposition+=1;

    oA=filename.substr(aposition,replacelength);
       //       G4cout << oA << G4endl;

    otarget=filename.substr(targetposition,3);
       //       G4cout << otarget << G4endl;

       
       // File names comparisons
       
    if ( (process == oprocess) || 
         (process == control)     ) processok = 1;
    else                            processok = 0;
       
    if (target == otarget)           targetok = 1;
    else                             targetok = 0;
       
    if (Z == oZ)                          Zok = 1;
    else
    {
      if (Z.substr(0,1) == control)       Zok = 1;
      else
      {
	oZvalue = atof(oZ);
	temp    = Z.substr(1,Z.length());
	Zvalue  = atof(temp);

	if (      (Z.substr(0,1) == ">" ) && (oZvalue > Zvalue) )  Zok = 1; 		 
	else if ( (Z.substr(0,1) == "<" ) && (oZvalue < Zvalue) )  Zok = 1; 		 
	else                                                       Zok = 0; 	 
      }
    }  
    if (A == oA)                  Aok = 1;
    else
    {
      if (A.substr(0,1)==control) Aok = 1;
      else
      {
	oAvalue = atof(oA);
	temp    = A.substr(1,A.length());
	Avalue  = atof(temp);

	if      ( (A.substr(0,1) == ">") && (oAvalue > Avalue) )  Aok = 1; 
		     
	else if ( (A.substr(0,1) == "<") && (oAvalue < Avalue) )  Aok = 1; 
		     
        else                                                      Aok = 0;      
      }   
    }   
    sumcheck = processok + targetok + Zok + Aok;
	   
    if ( sumcheck == 4 ) 
    {
      G4cout << "File found in" << G4endl;
      pathfilename = pathfile + "/" + originalfile;

      G4cout  << pathfilename << G4endl;
      fileout << pathfilename << std::endl;
      // fileout << "Check of fileout in StripFile method" << std::endl;
    }
    else
    {
       pathfilename="File not found";
       // fileout.close();
    }	   	   
  }
  fileout.close(); 
  errn = chdir(dirposition);
}


