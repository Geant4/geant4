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

#include "G4ReadHadDiffMomentumCXSC.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ElementVector.hh"

#include "g4std/fstream"
#include "g4std/strstream"
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


G4ReadHadDiffMomentumCXSC::G4ReadHadDiffMomentumCXSC( G4String&  primary,
                                              G4Isotope* isotope,
                                              G4String&  secondary   ):
  G4HadDataReading()
{
  ReadDiffMomentumCXSC( primary, isotope, secondary);
}

G4ReadHadDiffMomentumCXSC::G4ReadHadDiffMomentumCXSC( G4String&  primary,
                                              G4Element* element,
                                              G4String&  secondary   ):
  G4HadDataReading()
{
  ReadDiffMomentumCXSC( primary, element, secondary);
}

G4ReadHadDiffMomentumCXSC::G4ReadHadDiffMomentumCXSC( G4String&  primary,
                                              G4Material* material,
                                              G4String&  secondary   ):
  G4HadDataReading()
{
  ReadDiffMomentumCXSC( primary, material, secondary);
}



G4ReadHadDiffMomentumCXSC::~G4ReadHadDiffMomentumCXSC()
{ 

}



///////////////////////////////////////////////////////////////
//
//

void G4ReadHadDiffMomentumCXSC::ReadDiffMomentumCXSC( G4String&  primary,
                                         G4Isotope* isotope,
                                         G4String&  secondary   )
{
  // Build the complete string identifying the file with the data set
  
  G4int Z = isotope->GetZ();
  G4int N = isotope->GetN();
    
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);

  if( Z != 0 ) ost << "momZ"<< Z <<"N"<<N<< ".dat";  
  else         ost << "momN"<<N<<".dat"; 
  
  G4String name(nameChar);
  
  char* path = getenv("G4HADATASET");
  if ( !path )
  { 
    G4String excep = "G4ReadHadDiffMomentumCXSC - ";
    excep += "G4HADATASET environment variable not set";
    G4Exception(excep);
  }
  G4String pathString(path);
  G4String dirFile = pathString + "/" + primary + "/DxscDmomc/" + 
                                        secondary + "/" + name;
  G4bool momC = true;
  LoadDifferentialXSC(dirFile,momC);

}

///////////////////////////////////////////////////////////////
//
//

void G4ReadHadDiffMomentumCXSC::ReadDiffMomentumCXSC( G4String&  primary,
                                         G4Element* element,
                                         G4String&  secondary   )
{
  G4cout<<"G4ReadHadDiffMomentumCXSC::ReadDiffMomentumCXSC(p1,element,p2)"<<G4endl;
  
  G4int Z    = (G4int)element->GetZ();
  G4double A = element->GetN();     // A should be in xxx.xxx format 

  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  ost.setf(std::ios::fixed);
  ost.precision(3);  

  if( Z != 0 ) ost << "momZ"<< Z <<"A"<<A ;  
  else ;
  
  G4String name(nameChar);
  
  char* path = getenv("G4HADATASET");
  if ( !path )
  { 
    G4String excep = "G4ReadHadDiffMomentumCXSC - ";
    excep += "G4HADATASET environment variable not set";
    G4Exception(excep);
  }
  G4String pathString(path);
  G4String dirFile = pathString + "/" + primary + "/DxscDmomc/" + 
                                        secondary + "/" + name;
  G4bool momC = true;
  LoadDifferentialXSC(dirFile,momC);

}

///////////////////////////////////////////////////////////////
//
//

void G4ReadHadDiffMomentumCXSC::ReadDiffMomentumCXSC( G4String&  primary,
                                         G4Material* material,
                                         G4String&  secondary   )
{
  G4cout<<"G4ReadHadDiffMomentumCXSC::ReadDiffMomentumCXSC(p1,mat,p2)"<<G4endl;
  
  G4double Z = 0.;
  G4double A = 0.;

  for(size_t i = 0 ; i < material->GetNumberOfElements() ; ++i)
  {
    Z += material->GetElement(i)->GetZ()*
         material->GetFractionVector()[i];
    A += material->GetElement(i)->GetN()*
         material->GetFractionVector()[i];
  }
  // Z, A should be in xxx.xxx format  

  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  ost.setf(std::ios::fixed);
  ost.precision(3);  

  if( Z != 0 ) ost << "momZ"<< Z <<"A"<<A<<"mat" ;  
  else ;
  
  G4String name(nameChar);
  
  char* path = getenv("G4HADATASET");
  if ( !path )
  { 
    G4String excep = "G4ReadHadDiffMomentumCXSC - ";
    excep += "G4HADATASET environment variable not set";
    G4Exception(excep);
  }
  G4String pathString(path);
  G4String dirFile = pathString + "/" + primary + "/DxscDmomc/" + 
                                        secondary + "/" + name;
  G4bool momC = true;
  LoadDifferentialXSC(dirFile,momC);

}



//
//
/////////////////////////////////////////////////////////////










