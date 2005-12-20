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
// $Id: G4DNAAlphaPlusExcitation.cc,v 1.2 2005-12-20 13:52:19 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAAlphaPlusExcitation.hh"
#include "G4DNAGenericIonsManager.hh"

                                        G4DNAAlphaPlusExcitationEnergyLimitsPolicy :: G4DNAAlphaPlusExcitationEnergyLimitsPolicy()
:
 lowEnergyLimit(1*keV),
 zeroBelowLowEnergyLimit(false),
 highEnergyLimit(10*MeV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAAlphaPlusExcitationIncomingParticlePolicy :: G4DNAAlphaPlusExcitationIncomingParticlePolicy()
:
 kineticEnergyCorrection(0.9382723/3.727417)
{
 slaterEffectiveCharge[0]=2.0;
 slaterEffectiveCharge[1]=1.15;
 slaterEffectiveCharge[2]=1.15;
 sCoefficient[0]=0.7;
 sCoefficient[1]=0.15;
 sCoefficient[2]=0.15;
}

const G4ParticleDefinition *            G4DNAAlphaPlusExcitationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();
 
 return instance->GetIon("alpha+");
}
