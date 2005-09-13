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
// $Id: G4DNAProtonRuddIonization.cc,v 1.1 2005-09-13 08:59:45 zfrancis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAProtonRuddIonization.hh"
#include "G4Proton.hh"

                                        G4DNAProtonRuddIonizationEnergyLimitsPolicy :: G4DNAProtonRuddIonizationEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(false),
 highEnergyLimit(1000*keV),
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

