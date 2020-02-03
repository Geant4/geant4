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
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadeBuilder
//
// Author:      V.Ivanchenko 09.04.2006
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4IonBinaryCascadePhysics_h
#define G4IonBinaryCascadePhysics_h 1

#include "globals.hh"
#include "G4IonPhysics.hh"

class G4HadronInelasticProcess;
class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4VComponentCrossSection;
class G4FTFBuilder;
class G4BinaryLightIonReaction;

class G4IonBinaryCascadePhysics : public G4IonPhysics
{
public:
  G4IonBinaryCascadePhysics(G4int ver = 0);
  G4IonBinaryCascadePhysics(const G4String& name, G4int ver = 0);
  virtual ~G4IonBinaryCascadePhysics();

};


#endif

