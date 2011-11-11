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
// D. Wright, 23-Dec-2006 Added isotope dependence
// G.Folger, 25-Nov-2009: extend to 100TeV, using a constant above 20GeV
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element
//

#ifndef G4NeutronInelasticCrossSection_h
#define G4NeutronInelasticCrossSection_h

// Class Description
// Cross sections for neutron-nucleus scattering from 14 MeV up to 20 GeV,
// getting the low energy threshold behavior right.
// H.P. Wellisch (TRIUMF), M. Laidlaw (British Columbia U.). 1996. 
// Class Description - End

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"

class G4NeutronInelasticCrossSection : public G4VCrossSectionDataSet
{
public:

  G4NeutronInelasticCrossSection();
  ~G4NeutronInelasticCrossSection();
    
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle* aPart, 
			     G4int Z, const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  G4double GetCrossSection(G4double kineticEnergy, G4int Z, G4int A);

  virtual void CrossSectionDescription(std::ostream&) const;

private:

G4double minEnergy;
G4double maxEnergy;

};

#endif
