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
// $Id: G4DNAHydrogenRuddIonization.cc,v 1.2 2006-06-28 13:58:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAHydrogenRuddIonization.hh"
#include "G4DNAGenericIonsManager.hh"

                                        G4DNAHydrogenRuddIonizationEnergyLimitsPolicy :: G4DNAHydrogenRuddIonizationEnergyLimitsPolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(100*MeV),
 zeroAboveHighEnergyLimit(false)
{
}

                                        G4DNAHydrogenRuddIonizationIncomingParticlePolicy :: G4DNAHydrogenRuddIonizationIncomingParticlePolicy()

{
}

                                        G4DNAHydrogenRuddDataFilePolicy :: G4DNAHydrogenRuddDataFilePolicy()
:
 lowEnergyLimit(100*eV),
 zeroBelowLowEnergyLimit(true),
 highEnergyLimit(100*MeV),
 zeroAboveHighEnergyLimit(false),
 dataFileEnergyUnit(eV),
 dataFileCrossSectionUnit(m*m),
 dataFileName("RuddHydrogenIonizationCrossSection")
{
}

const G4ParticleDefinition *            G4DNAHydrogenRuddIonizationIncomingParticlePolicy :: IncomingParticleDefinition(void) const
{
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();

 return instance->GetIon("hydrogen");
}

