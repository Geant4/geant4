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
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 17-June-02
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4ChipsKaonZeroElasticXS for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
//
// ****************************************************************************************
// Short description: Cross-sections extracted (by W.Pokorski) from the CHIPS package for 
// K(zero)-nuclear  interactions. Original author: M. Kossov
// -------------------------------------------------------------------------------------
//

#include "G4ChipsKaonZeroElasticXS.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4KaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4AntiKaonZero.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsKaonZeroElasticXS);

// Initialization of the
G4double* G4ChipsKaonZeroElasticXS::lastLEN=0;// Pointer to the lastArray of LowEn CS
G4double* G4ChipsKaonZeroElasticXS::lastHEN=0;// Pointer to the lastArray of HighEnCS
G4int     G4ChipsKaonZeroElasticXS::lastN=0;  // The last N of calculated nucleus
G4int     G4ChipsKaonZeroElasticXS::lastZ=0;  // The last Z of calculated nucleus
G4double  G4ChipsKaonZeroElasticXS::lastP=0.; // Last used in cross section Momentum
G4double  G4ChipsKaonZeroElasticXS::lastTH=0.;// Last threshold momentum
G4double  G4ChipsKaonZeroElasticXS::lastCS=0.;// Last value of the Cross Section
G4int     G4ChipsKaonZeroElasticXS::lastI=0;  // The last position in the DAMDB


G4ChipsKaonZeroElasticXS::G4ChipsKaonZeroElasticXS():G4VCrossSectionDataSet(Default_Name())
{
  theKMinusCS = (G4ChipsKaonMinusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusElasticXS::Default_Name());;
  theKPlusCS = (G4ChipsKaonPlusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusElasticXS::Default_Name());
}


G4ChipsKaonZeroElasticXS::~G4ChipsKaonZeroElasticXS()
{
}

G4bool G4ChipsKaonZeroElasticXS::IsIsoApplicable(const G4DynamicParticle* Pt, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  G4ParticleDefinition* particle = Pt->GetDefinition();
  
  if(particle == G4KaonZeroLong::KaonZeroLong()   ||
     particle == G4KaonZeroShort::KaonZeroShort() ||
     particle == G4KaonZero::KaonZero()           ||
     particle == G4AntiKaonZero::AntiKaonZero()   )
  {
    return true;
  }
  return false;
}



// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)

G4double G4ChipsKaonZeroElasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
							const G4Isotope*,
							const G4Element*,
							const G4Material*)
{
  G4double pMom=Pt->GetTotalMomentum();
  G4int N = A - tgZ;
  
  return GetChipsCrossSection(pMom, tgZ, N, 311);
}

G4double G4ChipsKaonZeroElasticXS::GetChipsCrossSection(G4double mom, G4int Z, G4int N, G4int pdg)
{
  return (theKMinusCS->GetChipsCrossSection(mom,Z,N,pdg)
	  +theKPlusCS->GetChipsCrossSection(mom,Z,N,pdg))/2;
}


