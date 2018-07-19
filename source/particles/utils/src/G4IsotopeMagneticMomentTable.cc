//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4IsotopeMagneticMomentTable.cc
//
// Date:                16/03/07
// Author:              H.Kurashige
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
////////////////////////////////////////////////////////////////////////////////// IsomerLevel is added                         30 Apr. 2013  H.Kurashige

//
#include "G4IsotopeMagneticMomentTable.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include <fstream>
#include <sstream>

const G4double G4IsotopeMagneticMomentTable::levelTolerance = 2.0*keV;
// 0.1% torelance for excitation energy
  
const G4double G4IsotopeMagneticMomentTable::nuclearMagneton = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
// Nuclear Magneton
 
///////////////////////////////////////////////////////////////////////////////
G4IsotopeMagneticMomentTable::G4IsotopeMagneticMomentTable()
  :G4VIsotopeTable("MagneticMoment")
{
  if ( !getenv("G4IONMAGNETICMOMENT")) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IsotopeMagneticMomentTable::G4IsotopeMagneticMomentTable():  " 
	     <<  "Please setenv G4IONMAGNETICMOMENT for the magnetic moment data." 
	     << G4endl;
      G4Exception( "G4IsotopeMagneticMomentTable",
		   "File Not Found",
		   JustWarning, 
		   "Please setenv G4IONMAGNETICMOMENT");
    }
#endif
    G4Exception( "G4IsotopeMagneticMomentTable",
		 "File Not Found",
		 JustWarning, 
		 "Please setenv G4IONMAGNETICMOMENT");
    return;
  }
  
  G4String file = getenv("G4IONMAGNETICMOMENT");
  std::ifstream DataFile(file);

  if (!DataFile ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IsotopeMagneticMomentTable::G4IsotopeMagneticMomentTable():  " 
	    << file << " is not found " << G4endl;
    }
#endif
    G4Exception( "G4IsotopeMagneticMomentTable",
		 "File Not Found",
		 JustWarning, 
		 "Can not open G4IONMAGNETICMOMENT file");
    return;
  }
  
  char inputChars[80]={' '};
   
  while ( !DataFile.eof() ) { // Loop checking, 09.08.2015, K.Kurashige
    DataFile.getline(inputChars, 80);
    G4String inputLine = inputChars;
    G4int ionA, ionZ, ionJ, isomer;
    G4double ionE, ionMu, ionLife;
    G4String ionName, ionLifeUnit;
 
    if (inputChars[0] != '#' && inputLine.length() != 0) {
      std::istringstream tmpstream(inputLine);
      tmpstream >> ionZ >>  ionName >> ionA 
		>> isomer >>  ionE 
		>> ionLife >> ionLifeUnit
		>> ionJ >> ionMu;
          
      G4IsotopeProperty* fProperty = new G4IsotopeProperty();   
      // Set Isotope Property
      fProperty->SetAtomicNumber(ionZ);
      fProperty->SetAtomicMass(ionA);
      fProperty->SetIsomerLevel(isomer);
      fProperty->SetEnergy(ionE * MeV);
      fProperty->SetiSpin(ionJ);
      fProperty->SetMagneticMoment(ionMu*nuclearMagneton);
      
      fIsotopeList.push_back(fProperty);

      //if (GetVerboseLevel()>2) {
      // fProperty->DumpInfo();
      //}
	
    }
  }

  DataFile.close();
}

///////////////////////////////////////////////////////////////////////////////
G4IsotopeMagneticMomentTable::~G4IsotopeMagneticMomentTable()
{
  for (size_t i = 0 ; i< fIsotopeList.size(); i++) {
    delete fIsotopeList[i];
  }
  fIsotopeList.clear();
}

///////////////////////////////////////////////////////////////////////////////
G4IsotopeMagneticMomentTable::G4IsotopeMagneticMomentTable(const  G4IsotopeMagneticMomentTable & right)   
  :G4VIsotopeTable(right),
   fIsotopeList(0)
{
}

///////////////////////////////////////////////////////////////////////////////
G4IsotopeMagneticMomentTable & G4IsotopeMagneticMomentTable::operator= (const  G4IsotopeMagneticMomentTable &)
{
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
G4bool G4IsotopeMagneticMomentTable::FindIsotope(G4IsotopeProperty* pP)
{
  for (size_t i = 0 ; i< fIsotopeList.size(); ++i) {
    G4IsotopeProperty*  fP = fIsotopeList[i];
    
    // check Z
    if ( fP->GetAtomicNumber() > pP->GetAtomicNumber()) {
      // Not Found
      break;
    }
    if ( fP->GetAtomicNumber() < pP->GetAtomicNumber()) {
      // next
      continue;
    }
    
    // check A
    if ( fP->GetAtomicMass() != pP->GetAtomicMass()) {
      // next
      continue;
    }
    
    //check isomerLevel
    if (fP->GetIsomerLevel() != pP->GetIsomerLevel()) {
      // next
      continue;     
    }

    //check E
    if (std::fabs(fP->GetEnergy() - pP->GetEnergy()) < levelTolerance) {
      // Found
      return true;     
    }
    
  }
  return false;
}
///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* 
 G4IsotopeMagneticMomentTable::GetIsotope(G4int Z, G4int A, G4double E,
                                          G4Ions::G4FloatLevelBase /*flb*/)
{
  G4IsotopeProperty* fProperty = 0;
  for (size_t i = 0 ; i< fIsotopeList.size(); ++i) {
    G4IsotopeProperty*  fP = fIsotopeList[i];
 
     // check Z
    if ( fP->GetAtomicNumber() > Z) {
      // Not Found
      break;
    }
    if ( fP->GetAtomicNumber() < Z) {
      // next
      continue;
    }
    
    // check A
    if ( fP->GetAtomicMass() != A ) {
      // next
      continue;
    }
    
    //check E
    if (std::fabs(fP->GetEnergy() - E) < levelTolerance) {
      // Found
      fProperty = fP;
      // fP->DumpInfo();
      break;     
    }
    
  }

  return fProperty;
  
}

///////////////////////////////////////////////////////////////////////
G4IsotopeProperty* 
 G4IsotopeMagneticMomentTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl)
{
  G4IsotopeProperty* fProperty = 0;
  for (size_t i = 0 ; i< fIsotopeList.size(); ++i) {
    G4IsotopeProperty*  fP = fIsotopeList[i];
 
     // check Z
    if ( fP->GetAtomicNumber() > Z) {
      // Not Found
      break;
    }
    if ( fP->GetAtomicNumber() < Z) {
      // next
      continue;
    }
    // check A
    if ( fP->GetAtomicMass() != A ) {
      // next
      continue;
    }
    
    
    //check isomerLevel
    if (fP->GetIsomerLevel() == lvl) {
      // Found
      fProperty = fP;
      //fP->DumpInfo();
      break;     
    }
    
  }

  return fProperty;
  
}
