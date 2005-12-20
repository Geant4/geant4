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
// $Id: G4DNAProtonChargeDecrease.cc,v 1.2 2005-12-20 13:52:19 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAProtonChargeDecrease.hh"
#include "G4Proton.hh"

                                        G4DNAProtonChargeDecreaseEnergyLimitsPolicy :: G4DNAProtonChargeDecreaseEnergyLimitsPolicy()
:
 lowEnergyLimit(1.*keV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10.*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAProtonChargeDecreaseIncomingParticlePolicy :: G4DNAProtonChargeDecreaseIncomingParticlePolicy()
{
 f0[0]=1.;
 a0[0]=-0.180;
 a1[0]=-3.600;
 b0[0]=-18.22;
 b1[0]=-1.997;
 c0[0]=0.215;
 d0[0]=3.550;
 x0[0]=3.450;
 x1[0]=5.251;
}

const G4ParticleDefinition *            G4DNAProtonChargeDecreaseIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Proton::Proton();
}

G4int                                	G4DNAProtonChargeDecreaseIncomingParticlePolicy :: NumberOfPartialCrossSections(void) const
{
 return 1;
}
