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
// $Id: G4DNAProtonRuddIonization.cc,v 1.4 2006/06/29 19:39:42 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $

#include "G4DNAProtonRuddIonization.hh"
#include "G4Proton.hh"

                                        G4DNAProtonRuddIonizationEnergyLimitsPolicy :: G4DNAProtonRuddIonizationEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(500*keV),
 zeroAboveHighEnergyLimit(true)
{
}

                                        G4DNAProtonRuddIonizationIncomingParticlePolicy :: G4DNAProtonRuddIonizationIncomingParticlePolicy()

{
}

                                        G4DNAProtonRuddDataFilePolicy :: G4DNAProtonRuddDataFilePolicy()
:
 lowEnergyLimit(100*eV), 
 zeroBelowLowEnergyLimit(false),
 highEnergyLimit(1000*keV),
 zeroAboveHighEnergyLimit(true),
 dataFileEnergyUnit(eV),
 dataFileCrossSectionUnit(m*m),
 dataFileName("RuddProtonIonizationCrossSection")
{
}

const G4ParticleDefinition *            G4DNAProtonRuddIonizationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 return G4Proton::Proton();
}

