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
// $Id: G4NeutronHPArbitaryTab.hh,v 1.12 2007-06-06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPArbitaryTab_h
#define G4NeutronHPArbitaryTab_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4VNeutronHPEDis.hh"
#include "G4InterpolationManager.hh"

// we will need a List of these .... one per term.

class G4NeutronHPArbitaryTab : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPArbitaryTab()
  {
   theDistFunc = 0;
  }
  ~G4NeutronHPArbitaryTab()
  {
   if(theDistFunc!=0) delete [] theDistFunc;
  }
  
  inline void Init(std::ifstream & theData)
  {
    G4int i;
    theFractionalProb.Init(theData, eV);
    theData >> nDistFunc; // = number of incoming n energy points
    theDistFunc = new G4NeutronHPVector [nDistFunc];
    theManager.Init(theData);
    G4double currentEnergy;
    for(i=0; i<nDistFunc; i++)
    {
      theData >> currentEnergy;
      theDistFunc[i].SetLabel(currentEnergy*eV);
      theDistFunc[i].Init(theData, eV);
      theDistFunc[i].ThinOut(0.02); // @@@ optimization to be finished.
    }
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy) ;
  
  private:
  
  G4NeutronHPVector theFractionalProb;
  G4int nDistFunc;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4NeutronHPVector * theDistFunc; // one per incoming energy
  G4NeutronHPVector theBuffer;
  
};

#endif
