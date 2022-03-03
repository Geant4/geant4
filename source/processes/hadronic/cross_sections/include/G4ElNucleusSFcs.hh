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
// GEANT4 physics class: G4ElNucleusSFcs -- header file
//
//
// 30.09.21 V. Grichine eA inelastic cross section according to the structure function approach.
//                      The XS is calculated based on the corrected CHIPS XS.
//

#ifndef G4ElNucleusSFcs_h
#define G4ElNucleusSFcs_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include <vector>
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include <map>

class G4ElectroNuclearCrossSection;

class G4ElNucleusSFcs : public G4VCrossSectionDataSet
{
public:

  G4ElNucleusSFcs();
  virtual ~G4ElNucleusSFcs();
    
  static const char* Default_Name() {return "ElectronNucleusSFcs";}

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
                                     const G4Material*);
  
  virtual G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
                              const G4Isotope* iso = nullptr,
                              const G4Element* elm = nullptr,
			      const G4Material* mat = nullptr );
  G4double ThresholdEnergy();

  G4double GetRatio(G4int Z, G4int A);

private:

  static const G4double fZZ[19];
  static const G4double fRR[19];

  G4ElectroNuclearCrossSection* fCHIPScs;
};

#endif
