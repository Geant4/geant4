//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DNAHeliumChargeIncrease.cc,v 1.1 2005-12-20 13:41:32 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAHeliumChargeIncrease.hh"
#include "G4DNAGenericIonsManager.hh"

                                        G4DNAHeliumChargeIncreaseEnergyLimitsPolicy :: G4DNAHeliumChargeIncreaseEnergyLimitsPolicy()
:
 lowEnergyLimit(1.*keV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10.*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAHeliumChargeIncreaseIncomingParticlePolicy :: G4DNAHeliumChargeIncreaseIncomingParticlePolicy()
{
 f0[0]=1.;
 a0[0]=2.25;
 a1[0]=-0.75;
 b0[0]=-30.93;
 c0[0]=0.590;
 d0[0]=2.35;
 x0[0]=4.29;

 f0[1]=1.;
 a0[1]=2.25;
 a1[1]=-0.75;
 b0[1]=-32.61;
 c0[1]=0.435;
 d0[1]=2.70;
 x0[1]=4.45;

 // x1 and b1 will be calculated by G4DNADingfelderChargeChangeTotalCrossSectionPolicy<, > :: PartialCrossSection
 x1[0]=-1.;
 b1[0]=-1.;

 x1[1]=-1.;
 b1[1]=-1.;
}

const G4ParticleDefinition *            G4DNAHeliumChargeIncreaseIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();
 
 return instance->GetIon("helium");
}

G4int                                	G4DNAHeliumChargeIncreaseIncomingParticlePolicy :: NumberOfPartialCrossSections(void) const
{
 return 2;
}
