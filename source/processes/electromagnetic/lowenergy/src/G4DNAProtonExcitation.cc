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
// $Id: G4DNAProtonExcitation.cc,v 1.1 2005-07-20 10:01:54 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAProtonExcitation.hh"
#include "G4Proton.hh"

                                        G4DNAProtonExcitationEnergyLimitsPolicy :: G4DNAProtonExcitationEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(false),
 highEnergyLimit(100*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAProtonExcitationIncomingParticlePolicy :: G4DNAProtonExcitationIncomingParticlePolicy()
:
 kineticEnergyCorrection(1.)
{
 slaterEffectiveCharge[0]=0.;
 slaterEffectiveCharge[1]=0.;
 slaterEffectiveCharge[2]=0.;
 sCoefficient[0]=0.;
 sCoefficient[1]=0.;
 sCoefficient[2]=0.;
}

const G4ParticleDefinition *            G4DNAProtonExcitationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Proton::Proton();
}
