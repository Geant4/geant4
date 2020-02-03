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
#include "G4ParticleDefinition.hh"
#include "G4ConcreteNNToNNStar.hh"
#include "G4NNstarBuilder.hh"

G4ThreadLocal G4XNNstarTable *G4ConcreteNNToNNStar::theSigmaTable_G4MT_TLS_ = 0;

G4ConcreteNNToNNStar::G4ConcreteNNToNNStar(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary)  :
	G4ConcreteNNTwoBodyResonance(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
{
   if (!theSigmaTable_G4MT_TLS_) theSigmaTable_G4MT_TLS_ = new G4XNNstarTable;
   G4XNNstarTable &theSigmaTable = *theSigmaTable_G4MT_TLS_;
   establish_G4MT_TLS_G4ConcreteNNTwoBodyResonance(aPrimary,bPrimary,aSecondary,bSecondary,
		                                           G4NNstarBuilder(bSecondary->GetParticleName(),theSigmaTable));
}

G4ConcreteNNToNNStar::~G4ConcreteNNToNNStar()
{
	if (theSigmaTable_G4MT_TLS_) delete theSigmaTable_G4MT_TLS_;
	theSigmaTable_G4MT_TLS_=0;
}
