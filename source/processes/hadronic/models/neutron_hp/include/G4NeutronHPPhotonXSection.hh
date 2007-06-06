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
//
// $Id: G4NeutronHPPhotonXSection.hh,v 1.11 2007-06-06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPPhotonXSection_h
#define G4NeutronHPPhotonXSection_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4VNeutronVector.hh"

// we will need a List of these .... one per term.

class G4NeutronHPPhotonXSection 
{
  public:
  G4NeutronHPPhotonXSection()
  {
    theExclusive = 0;
    theExShell = 0;
    theExEnergy = 0;
    theExFlag = 0;
    theExDisFlag = 0;
  }
  ~G4NeutronHPPhotonXSection()
  {
    if(theExclusive!=0) delete [] theExclusive;
    if(theExShell != 0) delete [] theExShell;
    if(theExEnergy != 0) delete [] theExEnergy;
    if(theExFlag != 0) delete [] theExFlag;
    if(theExDisFlag != 0) delete [] theExDisFlag;
  }
  
  inline void Init(std::ifstream & aDataFile)
  {
    aDataFile  >> nChannels >> targetMass;
    if(nChannels!=1) 
    {
      aDataFile >> theIncEnergy>>theIncShell>>theIncFlag>>theIncDisFlag;
      theaDataFileInclusive.Init(aDataFile, eV);
    }
    theExclusive = new G4NeutronHPVector[nChannels];
    theExShell = new G4double[nChannels];
    theExEnergy = new G4double[nChannels];
    theExFlag = new G4int[nChannels];
    theExDisFlag = new G4int[nChannels];   
    for(G4int i=0; i<nChannels; i++)
    {
      aDataFile>>theExEnergy[i]>>theExShell[i]>>theExFlag[i]>>theExDisFlag[i];
      theExclusive[i].Init(aDataFile,eV);
    }
  }
  
  G4double Sample(G4double anEnergy)
  {
    return -1;
  }
  
  private:
   
  G4double targetMass;
  
  G4double theIncShell;
  G4double theIncEnergy;
  G4int theIncFlag
  G4int theIncDisFlag;
  G4NeutronHPVector theInclusive;
  
  G4int nChannels;
  G4double * theExShell;
  G4double * theExEnergy;
  G4int * theExFlag
  G4int * theExDisFlag;
  G4NeutronHPVector * theExclusive;
  
};

#endif
