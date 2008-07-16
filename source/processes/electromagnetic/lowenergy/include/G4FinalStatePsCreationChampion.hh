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
// $Id: G4FinalStatePsCreationChampion.hh,v 1.1 2008-07-16 19:01:07 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------


#ifndef G4FINALSTATEPSCREATIONCHAMPION_HH
#define G4FINALSTATEPSCREATIONCHAMPION_HH 1
 
#include "G4FinalStateProduct.hh"
#include "G4WaterIonisationStructure.hh"
#include "G4CrossSectionPsCreationChampionPartial.hh"
#include "G4DNAGenericIonsManager.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;

class G4FinalStatePsCreationChampion
{
public:
   
  G4FinalStatePsCreationChampion();
   
  ~G4FinalStatePsCreationChampion();
   
  const G4FinalStateProduct& GenerateFinalState(const G4Track& track, const G4Step& step);
   
private:

  G4String name;  
  G4double lowEnergyLimit;
  G4double highEnergyLimit;

  G4FinalStateProduct product;

  G4WaterIonisationStructure waterStructure;

  G4CrossSectionPsCreationChampionPartial cross;
};
#endif
