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

#include "globals.hh"
#include "G4ios.hh"
#include "G4XMesonBaryonElastic.hh"
#include "G4KineticTrack.hh"
#include "G4Gamma.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4XAqmElastic.hh"
#include "G4XPDGElastic.hh"

G4XMesonBaryonElastic::G4XMesonBaryonElastic()
{ 
  // As a first approximation the model is assumed to be valid over 
  // the entire energy range
  lowLimit = 0.;
  highLimit = DBL_MAX;
}


G4XMesonBaryonElastic::~G4XMesonBaryonElastic()
{ }


G4bool G4XMesonBaryonElastic::operator==(const G4XMesonBaryonElastic &right) const
{
  return (this == (G4XMesonBaryonElastic *) &right);
}


G4bool G4XMesonBaryonElastic::operator!=(const G4XMesonBaryonElastic &right) const
{
  return (this != (G4XMesonBaryonElastic *) &right);
}


G4double G4XMesonBaryonElastic::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4double sigma;

  // No gamma-baryon elastic scattering 
  const G4ParticleDefinition* defLight = FindLightParticle(trk1,trk2);
  if (defLight == G4Gamma::GammaDefinition())
    {
      sigma = 0;
    }
  else
    {

      G4LorentzVector p41 = trk1.Get4Momentum();
      G4ThreeVector p3 = trk1.GetPosition();
      G4ParticleDefinition* def = G4PionPlus::PionPlusDefinition();

      G4KineticTrack piTrk(def,
			   trk1.GetFormationTime(),
			   p3,
			   (G4LorentzVector&)p41);

      G4LorentzVector p42 = trk2.Get4Momentum();    
      G4KineticTrack pTrk(((G4ParticleDefinition*)G4Proton::ProtonDefinition()),
			  trk2.GetFormationTime(),
			  (G4ThreeVector)trk2.GetPosition(),
			  (G4LorentzVector&)p42);
      
      G4XAqmElastic aqm;
      G4double xAqmDummy = aqm.CrossSection(piTrk,pTrk);
      G4double xAqm = aqm.CrossSection(trk1,trk2);
      G4double factor = 1.;
      if (xAqmDummy != 0.0)
	{
	  factor = xAqm / xAqmDummy;
	}
      G4XPDGElastic pdg;
      
      sigma = pdg.CrossSection(piTrk,pTrk);
      sigma = sigma * factor;
    }
  
  return sigma;
}


G4String G4XMesonBaryonElastic::Name() const
{
  G4String name("MesonBaryonElasticCrossSection");
  return name;
}



G4bool G4XMesonBaryonElastic::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,lowLimit,highLimit);

  return answer;
}
