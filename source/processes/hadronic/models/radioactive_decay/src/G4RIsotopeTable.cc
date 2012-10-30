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
// MODULE:              G4RIsotopeTable.cc
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 14 April 2000, F Lei, DERA UK
// 0.b.4 release. Minor changes to 
//            1) levelTolerance = 2.0 keV
//            2) changes to verbose control
//
// 18,July 2001 F.Lei
//  tidy up the print out at run level
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4DecayTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IsotopeProperty.hh"
#include "G4RIsotopeTable.hh"

#include "G4HadronicException.hh"
#include "G4NuclearLevelStore.hh"

/*
#include "G4RadioactiveDecayMode.hh"
#include "G4ITDecayChannel.hh"
#include "G4BetaMinusDecayChannel.hh"
#include "G4BetaPlusDecayChannel.hh"
#include "G4KshellECDecayChannel.hh"
#include "G4LshellECDecayChannel.hh"
#include "G4AlphaDecayChannel.hh"
*/
#include "G4ios.hh"
#include "globals.hh"
#include <iomanip>
#include <fstream>
#include <sstream>

const G4double G4RIsotopeTable::levelTolerance = 2.0*keV;

///////////////////////////////////////////////////////////////////////////////
//
G4RIsotopeTable::G4RIsotopeTable()
{//
 //Reset the list of user define data file
 //
 theUserRadioactiveDataFiles.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
G4RIsotopeTable::~G4RIsotopeTable()
{
  fIsotopeList.clear();
  fIsotopeNameList.clear(); 
}
///////////////////////////////////////////////////////////////////////////////
//
G4int G4RIsotopeTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool G4RIsotopeTable::FindIsotope(G4IsotopeProperty* )
{
  // do nothing, it is here just for the compiler
  // it is required by the base class
  return true;
}
///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* G4RIsotopeTable::GetIsotope(G4int Z, G4int A, G4double E)
{
  G4String fname = GetIsotopeName(Z, A, E);  
  G4int j = -1;
  for (G4int i = 0 ; i< Entries(); i++) {
    if(fIsotopeNameList[i] == fname) j = i;}
  if (j >=0) {
    if (GetVerboseLevel()>1) {
    G4cout <<"G4RIsotopeTable::GetIsotope No. : ";
    G4cout <<j<<G4endl;   
    }
    return  GetIsotope(j);}
  // isotope property data has been loaded already and just return the pointer
  else{
    G4double meanlife = GetMeanLifeTime(Z, A, E);
    // E is pass as a refence hence on entry E is supplied by the user and it 
    // could be slightly different from the returned value which is the one 
    // defined in the database.
    // this call is to ensure the code uses a consistane E value through out.
    //
    
    G4IsotopeProperty* fProperty = new G4IsotopeProperty();   
    // Set Isotope Property
    fProperty->SetLifeTime(meanlife);
    fProperty->SetAtomicNumber(Z);
    fProperty->SetAtomicMass(A);
    // Notic that the value of E may have been changed
    fProperty->SetEnergy(E);
    // The spin is not being used in the current implementation
    fProperty->SetiSpin(0);
    // the decaytable will be loaded later in G4RadioactiveDecay when it is needed
    fProperty->SetDecayTable(0);
    
    fIsotopeList.push_back(fProperty);
    fname = GetIsotopeName(Z, A, E);
    fIsotopeNameList.push_back(fname);
    if (GetVerboseLevel()>1) {
      G4cout <<"G4RIsotopeTable::GetIsotope create: ";
      G4cout <<fname <<G4endl;  
    }
    return fProperty;

  }
}
///////////////////////////////////////////////////////////////////////////////
//
G4String G4RIsotopeTable::GetIsotopeName(G4int Z, G4int A, G4double E)  
{
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os <<"A"<< A << "Z" << Z <<'[' << std::setprecision(1) << E/keV << ']';
  G4String name = os.str();
  if (GetVerboseLevel()>1) {
    G4cout <<"G4RIsotopeTable::GetIsotope Name: ";
    G4cout <<name <<G4endl;   
  }
  return name;
}


G4double G4RIsotopeTable::GetMeanLifeTime(G4int Z, G4int A, G4double& aE)
{

  G4double lifetime = -1.0;


  //Check if data have been provided by the user
  G4String file= theUserRadioactiveDataFiles[1000*A+Z];
  if (file ==""){
	if (!getenv("G4RADIOACTIVEDATA")) {
		G4cout << "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files." << G4endl;
		throw G4HadronicException(__FILE__, __LINE__,
			      "Please setenv G4RADIOACTIVEDATA to point to the radioactive decay data files.");
	}
	G4String dirName = getenv("G4RADIOACTIVEDATA");

	std::ostringstream os;
	os <<dirName <<"/z" <<Z <<".a" <<A ;
	file = os.str();
  }
  std::ifstream DecaySchemeFile(file);

  G4bool found_in_raddecay_data(false);
  if (!DecaySchemeFile) {
    if (GetVerboseLevel()>1) {
      G4cout <<"G4RIsotopeTable::GetMeanLife() : "
	     <<"cannot find ion radioactive decay file: " 
	     <<file <<G4endl;
    }
  } else {
    char inputChars[100]={' '};
    G4String inputLine;
    G4String recordType("");
    G4double a(0.0);
    G4double b(0.0);

    while (!found_in_raddecay_data && !DecaySchemeFile.getline(inputChars, 100).eof()) {
      inputLine = inputChars;
      inputLine = inputLine.strip(1);

      if (inputChars[0] != '#' && inputLine.length() != 0) {
        std::istringstream tmpstream(inputLine);
        tmpstream >> recordType >> a >> b;
        if (recordType == "P") {
          if (std::abs(a*keV-aE) < levelTolerance) {
            found_in_raddecay_data    = true;
            lifetime = b/0.693147*s ;
          }
        }
      }
    }
    DecaySchemeFile.close();
  }

    if (!found_in_raddecay_data && aE) {
      G4double half_life=-1.;
      lifetime = 1.0E-20*s;


      //added by L.Desorgher If the life time is not found in  raddecay database
      // then it is deduced from photo-evaporation level
      const G4NuclearLevel* aLevel =
    		G4NuclearLevelStore::GetInstance()->GetManager(Z, A)
    									->NearestLevel(aE,levelTolerance);
      if (aLevel) {
    	  half_life = aLevel->HalfLife();
    	  lifetime = half_life/0.693147;
      }

      if (GetVerboseLevel()>1 && half_life<0) {
        G4cout << "G4RIsotopeTable::GetMeanLife() : ";
        G4cout << "cannot find ion of required excitation E = " << aE << G4endl;
        G4cout << "state in radioactive or photoevaporation data file " << G4endl;
        G4cout <<"The nucleus is assumed to be IT decayed with life = 1E-20 s" << G4endl;
        G4cout <<" -----------* THIS MAY CAUSE PROBLEM IN ITS DECAY-----------" << G4endl;
      }
    }

    if (!found_in_raddecay_data && !aE) {
      if (GetVerboseLevel()>1) {
        G4cout <<"G4RIsotopeTable::GetMeanLife() : ";
        G4cout <<"cannot find ion of required excitation E = " << aE << G4endl;
        G4cout <<"state in radioactive or photoevaporation data file" <<G4endl;
        G4cout <<"The nucleus is assumed to be stable" <<G4endl;
        lifetime = -1.0;
      }
    }

    if (GetVerboseLevel()>1) {
       G4cout <<"G4RIsotopeTable::GetMeanLifeTime: ";
       G4cout <<lifetime << " for " << GetIsotopeName(Z, A, aE) <<G4endl;
    }
  return lifetime;
}
////////////////////////////////////////////////////////////////////
//
void G4RIsotopeTable::AddUserDecayDataFile(G4int Z, G4int A,G4String filename)
{ if (Z<1 || A<2) {
	G4cout<<"Z and A not valid!"<<G4endl;
  }

  std::ifstream DecaySchemeFile(filename);
  if (DecaySchemeFile){
	G4int ID_ion=A*1000+Z;
	theUserRadioactiveDataFiles[ID_ion]=filename;
  }
  else {
	G4cout<<"The file "<<filename<<" does not exist!"<<G4endl;
  }
}

