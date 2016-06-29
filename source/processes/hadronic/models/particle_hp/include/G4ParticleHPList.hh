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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPList_h
#define G4ParticleHPList_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

class G4ParticleHPList
{
  public:
  
  G4ParticleHPList()
  {
    theData = new G4double[2]; 
    nPoints=2;
    nEntries=0;
    theLabel=0.0;
  }
  
  ~G4ParticleHPList()
  {
    delete [] theData;
  }
  
  inline void SetValue(G4int i, G4double y) 
  { 
    Check(i);
    theData[i]=y;
  }
  G4double GetValue(G4int i);
  
  inline G4int GetListLength() {return nEntries;}

  void Dump();
  
  void Init(std::istream & aDataFile, G4int nPar, G4double unit=1.);
  
  void Init(std::istream & aDataFile, G4double unit=1.);

  inline void SetLabel(G4double aLabel) { theLabel = aLabel; }
  
  inline G4double GetLabel() { return theLabel; }

  private:
  
  void Check(G4int i);
 
  G4double theLabel;

  G4double * theData;
  G4int nEntries;
  G4int nPoints;
};

#endif
