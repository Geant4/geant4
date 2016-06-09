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
// $Id: G4DNAProtonBornExcitation.cc,v 1.3 2006/06/29 19:39:34 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

#include "G4DNAProtonBornExcitation.hh"
#include "G4Proton.hh"

                                        G4DNAProtonBornExcitationEnergyLimitsPolicy :: G4DNAProtonBornExcitationEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(100*MeV),
 zeroAboveHighEnergyLimit(false)
{
}

                                        G4DNAProtonBornExcitationIncomingParticlePolicy :: G4DNAProtonBornExcitationIncomingParticlePolicy()

{
}

                                        G4DNAProtonBornExcitationDataFilePolicy :: G4DNAProtonBornExcitationDataFilePolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(100*MeV),
 zeroAboveHighEnergyLimit(false),
 dataFileEnergyUnit(eV),
 dataFileCrossSectionUnit(m*m),
 dataFileName("ProtonBornExcitationCrossSection")
{
}

const G4ParticleDefinition *            G4DNAProtonBornExcitationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Proton::Proton();
}

