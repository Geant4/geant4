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
// $Id: G4DNAProtonChargeDecrease.cc,v 1.1 2005-09-15 09:04:21 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAProtonChargeDecrease.hh"
#include "G4Proton.hh"

                                        G4DNAProtonChargeDecreaseEnergyLimitsPolicy :: G4DNAProtonChargeDecreaseEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(2*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAProtonChargeDecreaseIncomingParticlePolicy :: G4DNAProtonChargeDecreaseIncomingParticlePolicy()

{}

const G4ParticleDefinition *            G4DNAProtonChargeDecreaseIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Proton::Proton();
}
