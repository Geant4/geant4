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
// $Id: G4DiffElasticRatio.cc 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4DiffElasticRatio
//
// Author:  V.Grichine 07.10.2014
//
// Modifications:
//
// 16.03.15 V. Grichine safety against H ( A > 1 only ) 

#include "G4DiffElasticRatio.hh"
#include "G4ParticleDefinition.hh"
#include "G4DiffElasticRatio.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"

G4DiffElasticRatio::G4DiffElasticRatio(const G4String& nam, G4int verb) 
  : G4VCrossSectionRatio( nam, verb) 
{
  fGGXsc = new G4ComponentGGHadronNucleusXsc();
  fDDthreshold = 450.*CLHEP::MeV; // ~3 pi masses
}

G4DiffElasticRatio::~G4DiffElasticRatio()
{
  if(fGGXsc) delete fGGXsc;
}


G4double G4DiffElasticRatio::ComputeRatio(const G4ParticleDefinition* theParticleDefinition,
			G4double kinEnergy, 
			G4int Z, G4int A)
{
  G4double ratio = 0.;

  if( A > 1 && kinEnergy > fDDthreshold )
  {
    G4double ggElXsc = fGGXsc->GetElasticElementCrossSection(theParticleDefinition,kinEnergy,
                                                       Z,A);
    G4double ggsdXsc = fGGXsc->GetDiffractionGlauberGribovXsc();

    if( ggElXsc > 0.) ratio = ggsdXsc/ggElXsc;
    else              ratio = 0;
  }
  // G4cout<<theParticleDefinition->GetParticleName()<<"; "<<kinEnergy/CLHEP::GeV<<" GeV; r = "<<ratio<<G4endl;
  return ratio;
}
