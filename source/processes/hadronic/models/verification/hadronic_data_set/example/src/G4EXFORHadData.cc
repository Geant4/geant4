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

#include "G4EXFORHadData.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ElementVector.hh"

#include <fstream>
#include <strstream>
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//sim
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "G4HadFileFinder.hh"
//end sim

G4String G4EXFORHadData::EXFORenergy = "MeV";
G4String G4EXFORHadData::EXFORangle = "degree";
G4String G4EXFORHadData::EXFORddXsc = "barn/sr/MeV";

G4EXFORHadData::G4EXFORHadData( G4String&  inputFn, G4HadFileSpec& filetowrite):
  G4VHadDataWriting()
{
  SetInputFileName(inputFn);
  SetEnergyUnit(EXFORenergy);
  SetAngleUnit(EXFORangle);
  SetDdXscUnit(EXFORddXsc);

  WriteDataFile(filetowrite);
}

G4EXFORHadData::~G4EXFORHadData()
{ 

}



///////////////////////////////////////////////////////////////
//
//

void G4EXFORHadData::WriteDataFile( G4HadFileSpec& filetowrite) 
{

  G4String name(filetowrite.G4HDSFilename());
  G4String pathString(filetowrite.G4HDSFilepath());
  
  G4String dirFile = pathString+name;  

  char* doIhaverights = getenv("G4HADWRITTINGRIGHTS");

  int doI=strcmp(doIhaverights,"0");

  if ((!doIhaverights)||(doI==0))
    {
      G4String excep = "G4HADWRITTINGRIGHTS environment variable not set or 0. ";
      excep += "You don't have the permission to write in the G4HDS";
      G4Exception(excep);
      
    }


  G4cout << "Filling file in position " << dirFile << G4endl;

  FillDoubleDiffXSC(dirFile);

}











