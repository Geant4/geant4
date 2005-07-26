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
// $Id: G4DNAElectronExcitation.cc,v 1.1 2005-07-26 12:26:13 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAElectronExcitation.hh"
#include "G4Electron.hh"

                                        G4DNAElectronExcitationEnergyLimitsPolicy :: G4DNAElectronExcitationEnergyLimitsPolicy()
:
 lowEnergyLimit(7*eV),
 zeroBelowLowEnergyLimit(false),
 highEnergyLimit(10*keV),
 zeroAboveHighEnergyLimit(false)
{
}

                                        G4DNAElectronExcitationIncomingParticlePolicy :: G4DNAElectronExcitationIncomingParticlePolicy()

{}

const G4ParticleDefinition *            G4DNAElectronExcitationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Electron::Electron();
}
