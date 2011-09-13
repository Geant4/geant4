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
// by JPW, working, but to be cleaned up. @@@@
// G.Folger, 29-sept-2006: extend to 1TeV, using a constant above 20GeV
// D. Wright, 23-Dec-2006: added isotope dependence
// G.Folger, 25-Nov-2009: extend to 100TeV, using a constant above 20GeV
// V.Ivanchenko, 18-Aug-2011: migration to new design
//

#ifndef G4ProtonInelasticCrossSection_h
#define G4ProtonInelasticCrossSection_h

// Class Description
// Cross-sections for proton nuclear scattering up to 20 GeV, getting the low
// energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), D. Axen (British Columbia U.). 1996. 
// Published in Phys.Rev.C54:1329-1332,1996
// Implements corrected parameterization from  
// http://laws.lanl.gov/XCI/PEOPLE/rep/pdf/RPS96.pdf
// Class Description - End


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4NistManager;
class G4Proton;

class G4ProtonInelasticCrossSection : public G4VCrossSectionDataSet
{
public:

  G4ProtonInelasticCrossSection();

  ~G4ProtonInelasticCrossSection();

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle* aPart, 
			     G4int Z, const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  G4double GetProtonCrossSection(G4double kineticEnergy, G4int Z);

private:

  G4NistManager* nist;
  G4Proton* theProton;
  G4double thEnergy;

};

#endif
