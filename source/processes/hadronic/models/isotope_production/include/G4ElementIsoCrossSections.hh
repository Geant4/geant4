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
#ifndef G4ElementIsoCrossSections_h
#define G4ElementIsoCrossSections_h

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
