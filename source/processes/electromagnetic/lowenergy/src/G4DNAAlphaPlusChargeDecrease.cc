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
// $Id: G4DNAAlphaPlusChargeDecrease.cc,v 1.2 2006/06/29 19:38:50 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $

#include "G4DNAAlphaPlusChargeDecrease.hh"
#include "G4DNAGenericIonsManager.hh"

                                        G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy :: G4DNAAlphaPlusChargeDecreaseEnergyLimitsPolicy()
:
 lowEnergyLimit(1.*keV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10.*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy :: G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy()
{
 f0[0]=1.;
 a0[0]=0.65;
 a1[0]=-2.75;
 b0[0]=-21.81;
 c0[0]=0.232;
 d0[0]=2.95;
 x0[0]=3.53;

 // x1 and b1 will be calculated by G4DNADingfelderChargeChangeTotalCrossSectionPolicy<, > :: PartialCrossSection
 x1[0]=-1.;
 b1[0]=-1.;
}

const G4ParticleDefinition *            G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();
 
 return instance->GetIon("alpha+");
}

G4int                                	G4DNAAlphaPlusChargeDecreaseIncomingParticlePolicy :: NumberOfPartialCrossSections(void) const
{
 return 1;
}
