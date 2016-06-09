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
// $Id: G4DNAAlphaPlusChargeIncrease.cc,v 1.2 2006/06/29 19:38:52 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

#include "G4DNAAlphaPlusChargeIncrease.hh"
#include "G4DNAGenericIonsManager.hh"

                                        G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy :: G4DNAAlphaPlusChargeIncreaseEnergyLimitsPolicy()
:
 lowEnergyLimit(1.*keV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10.*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy :: G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy()
{
 f0[0]=1.;
 a0[0]=2.25;
 a1[0]=-0.75;
 b0[0]=-32.10;
 c0[0]=0.600;
 d0[0]=2.40;
 x0[0]=4.60;

 // x1 and b1 will be calculated by G4DNADingfelderChargeChangeTotalCrossSectionPolicy<, > :: PartialCrossSection
 x1[0]=-1.;
 b1[0]=-1.;
}

const G4ParticleDefinition *            G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();
 
 return instance->GetIon("alpha+");
}

G4int                                	G4DNAAlphaPlusChargeIncreaseIncomingParticlePolicy :: NumberOfPartialCrossSections(void) const
{
 return 1;
}
