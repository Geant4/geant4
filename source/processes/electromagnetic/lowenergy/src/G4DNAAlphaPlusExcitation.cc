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
// $Id: G4DNAAlphaPlusExcitation.cc,v 1.3 2006/06/29 19:38:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

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
