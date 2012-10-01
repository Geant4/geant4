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
#ifndef G4ElementIsoCrossSections_h
#define G4ElementIsoCrossSections_h

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4StableIsotopes.hh"
#include "G4IsoResult.hh"
#include "Randomize.hh"

template <class IsoIsoCrossSectionType>
class G4ElementIsoCrossSections
{

public:
  
  G4ElementIsoCrossSections();
  ~G4ElementIsoCrossSections();
  void Init(const G4Element * anElement);
  
  G4double GetCrossSection(G4double anEnergy);
  G4IsoResult * GetProductIsotope(G4double anEnergy);

private:
  
  IsoIsoCrossSectionType ** theData;
  G4int nIsotopes;
  G4StableIsotopes theStableOnes;
  G4double crossSectionBuffer;

};

#include "G4ElementIsoCrossSections.icc"

#endif
