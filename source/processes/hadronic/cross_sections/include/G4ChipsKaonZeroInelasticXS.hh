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
//
//
// GEANT4 physics class: G4ChipsKaonZeroInelasticXS -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//
// ****************************************************************************************
// Short description: Cross-sections extracted (by W.Pokorski) from the CHIPS package for 
// K(zero)-nuclear  interactions. Original author: M. Kossov
// -------------------------------------------------------------------------------------
//

#ifndef G4ChipsKaonZeroInelasticXS_h
#define G4ChipsKaonZeroInelasticXS_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "G4VCrossSectionDataSet.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"

class G4ChipsKaonZeroInelasticXS : public G4VCrossSectionDataSet
{


public:

  G4ChipsKaonZeroInelasticXS();

  ~G4ChipsKaonZeroInelasticXS();

  static const char* Default_Name() {return "ChipsKaonZeroInelasticXS";}

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual G4bool IsIsoApplicable(const G4DynamicParticle* Pt, G4int Z, G4int A,    
				 const G4Element* elm,
				 const G4Material* mat );

  // At present momentum (pMom) in MeV/c, CS in mb (@@ Units)
  virtual G4double GetIsoCrossSection(const G4DynamicParticle*, G4int tgZ, G4int A,  
				      const G4Isotope* iso = 0,
				      const G4Element* elm = 0,
				      const G4Material* mat = 0);

  G4double GetChipsCrossSection(G4double momentum, G4int Z, G4int N, G4int pdg);  

// Body
private:

  G4ChipsKaonMinusInelasticXS* theKMinusCS; // K- cross-section
  G4ChipsKaonPlusInelasticXS* theKPlusCS;  // K+ cross-section

  G4double* lastLEN; // Pointer to the last array of LowEnergy cross sections
  G4double* lastHEN; // Pointer to the last array of HighEnergy cross sections
  G4int     lastN;   // The last N of calculated nucleus
  G4int     lastZ;   // The last Z of calculated nucleus
  G4double  lastP;   // Last used in the cross section Momentum
  G4double  lastTH;  // Last value of the Momentum Threshold
  G4double  lastCS;  // Last value of the Cross Section
  G4int     lastI;   // The last position in the DAMDB

};

#endif
