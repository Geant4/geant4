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
#include <strstream>
#include <iostream>
#include "G4ios.hh"
#include "globals.hh"


#include "G4HadFileSpec.hh"


////////////////////////////////////////////////////////////////
//
//

G4HadFileSpec::G4HadFileSpec(G4String&  primary, G4Isotope* isotope, G4String&  secondary, G4String& process)
  :fprimary(primary),fisotope(isotope),fsecondary(secondary),fprocess(process)
{
  fmaterial=NULL;
  felement=NULL;
}


G4HadFileSpec::G4HadFileSpec(G4String&  primary, G4Element* element, G4String&  secondary, G4String& process)
  :fprimary(primary),felement(element),fsecondary(secondary),fprocess(process)
{
  fisotope=NULL;
  fmaterial=NULL;
}

G4HadFileSpec::G4HadFileSpec(G4String&  primary, G4Material* material, G4String&  secondary, G4String& process)
  :fprimary(primary),fmaterial(material),fsecondary(secondary),fprocess(process)
{
  fisotope=NULL;
  felement=NULL;
}



////////////////////////////////////////////////////////////////
//
//

G4String G4HadFileSpec::G4HDSFilename() 
{

  char nameChar[100] = {""};
  if (fisotope!=NULL)
    {
      G4int Zi = fisotope->GetZ();
      G4int Ni = fisotope->GetN();
  
      std::ostrstream ost(nameChar, 100, std::ios::out);

      if( Zi != 0 ) ost <<fprocess<< "Z"<< Zi <<".000N"<<Ni<<".000"<<"iso";  
      else         ost <<fprocess<< "N"<<Ni<<".000"<<"iso"; 
    }
  else if(felement!=NULL)
    {
      G4int Ze    = (G4int)felement->GetZ();
      G4double Ae = felement->GetN();     // A should be in xxx.xxx format 
      
      std::ostrstream ost(nameChar, 100, std::ios::out);
      ost.setf(std::ios::fixed);
      ost.precision(3);  
      
      if( Ze != 0 ) ost<<fprocess<< "Z"<< Ze <<".000A"<<Ae<<"ele";  
      else ;
    }
  else if(fmaterial!=NULL)
   {
     G4double Zm = 0.;
     G4double Am = 0.;
     
     for(size_t i = 0 ; i < fmaterial->GetNumberOfElements() ; ++i)
       {
	 Zm += fmaterial->GetElement(i)->GetZ()*
	   fmaterial->GetFractionVector()[i];
	 Am += fmaterial->GetElement(i)->GetN()*
	   fmaterial->GetFractionVector()[i];
       }
     // Z, A should be in xxx.xxx format  
     
     std::ostrstream ost(nameChar, 100, std::ios::out);
     ost.setf(std::ios::fixed);
     ost.precision(3);  

     if( Zm != 0 ) ost <<fprocess<< "Z"<< Zm <<"A"<<Am<<"mat" ;  
     else ;
   }
     

  G4String name(nameChar);

  return name;
}


G4String G4HadFileSpec::G4HDSFilepath()
{

  G4String fprocessname;



  char* path = getenv("G4HADATASET");
  if ( !path )
    { 
      G4String excep = "";
      excep = "G4HADATASET environment variable not set";
      G4Exception(excep);
    };
  
  if(fprocess=="dd")
    {
      fprocessname="doublediff";
    }
  
  else if(fprocess=="df")
    {
      fprocessname="diff";
    }
  
  else if(fprocess=="el")
    {
      fprocessname="elastic";
    }
  else if(fprocess=="in")
    {
      fprocessname="inelastic";
    }
  else if(fprocess=="mu")
    {
      fprocessname="multeplicity";
    }
  else if(fprocess=="to")
    {
      fprocessname="total";
    }
  else
    {
      G4String excep = "";
      G4String mm = "The process you are asking for is not in the G4HDS predifined list.";
      G4cout << mm << G4endl;
      excep+="Check the list in $G4HADATASET/G4HDSProcesses.txt.";
      G4Exception(excep);
    };
  


  G4String pathString(path);
  G4String dirFile = pathString + "/" + fprimary + "/" + fprocessname +"/" + 
    fsecondary + "/";

    
  if (chdir(dirFile)==-1)
    {
      G4String whatYouwant;
      
      G4cout << " The directory you are asking for doesn't exist " << G4endl;
      G4cout << " Do you want me to create it (y/n)?"<< G4endl;
      G4cin >> whatYouwant;

      char* doIhaverights = getenv("G4HADWRITTINGRIGHTS");
      
      int doI=strcmp(doIhaverights,"0");
      
      if (((!doIhaverights)||(doI==0))&&(whatYouwant=="y"))
	{ 
	  G4String message = "G4HADWRITTINGRIGHTS environment variable not set or 0";
	  G4String excep = "You have only reading permission in G4HDS";
	  G4cout << message << G4endl;
	  G4Exception(excep);      
	}
      else if((doI==1)&&(whatYouwant=="y"))
	{
	  G4String message3 = "Warning: you have writting permission in G4HDS";
	  G4cout << message3 << G4endl;   
	  G4String message4 = "I am creating the directory";
	  G4cout << message4 << G4endl; 

	  G4String dirname=pathString + "/" + fprimary;
	  int errornum=mkdir(dirname,0777);
	  dirname+="/"+fprocessname+"/";
	  errornum=mkdir(dirname,0777);
	  dirname+=fsecondary;
	  errornum=mkdir(dirname,0777);

	}
      else if(((!doIhaverights)||(doI==0))&&(whatYouwant=="n"))
	{
	  G4String excep = "So I stop here";
	  G4Exception(excep);
	}
     else if((doI==1)&&(whatYouwant=="n"))
	{
	  G4String excep = "So I stop here";
	  G4Exception(excep);
	} 
	
  }
  // int err = chdir(getenv("G4HADWORKING"));
  chdir(getenv("G4HADWORKING"));
  
  return dirFile;
}
