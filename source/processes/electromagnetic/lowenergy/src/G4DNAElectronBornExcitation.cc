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
// $Id: G4DNAElectronBornExcitation.cc,v 1.2 2006-06-28 13:58:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAElectronBornExcitation.hh"
#include "G4Electron.hh"

                                        G4DNAElectronBornExcitationEnergyLimitsPolicy::G4DNAElectronBornExcitationEnergyLimitsPolicy()
:
 lowEnergyLimit(7*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10*keV),
 zeroAboveHighEnergyLimit(false)
{
}

                                        G4DNAElectronBornExcitationIncomingParticlePolicy :: G4DNAElectronBornExcitationIncomingParticlePolicy()

{
}

                                        G4DNAElectronBornExcitationDataFilePolicy :: G4DNAElectronBornExcitationDataFilePolicy()
:
 lowEnergyLimit(7*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(10*keV),
 zeroAboveHighEnergyLimit(false),
 dataFileEnergyUnit(eV),
 dataFileCrossSectionUnit(m*m),
 dataFileName("ElectronBornExcitationCrossSection")
{
}

const G4ParticleDefinition *            G4DNAElectronBornExcitationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Electron::Electron();
}

