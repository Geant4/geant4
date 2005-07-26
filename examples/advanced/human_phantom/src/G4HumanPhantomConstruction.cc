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
#include "globals.hh"

#include "G4HumanPhantomConstruction.hh"
#include "G4HumanPhantomSD.hh"

#include "G4VBodyFactory.hh"
#include "G4MIRDBodyFactory.hh"
#include "G4ORNLBodyFactory.hh"

#include "G4FemaleBuilder.hh"
#include "G4MaleBuilder.hh"

#include "G4SDManager.hh"

//#include "G4VisAttributes.hh"   
//#include "G4Colour.hh"


G4HumanPhantomConstruction::G4HumanPhantomConstruction()
{
}

G4HumanPhantomConstruction::~G4HumanPhantomConstruction()
{
}


G4VPhysicalVolume* G4HumanPhantomConstruction::Construct()
{

//   G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
//   logicWorld  ->SetVisAttributes(BoxVisAtt);

  G4FemaleBuilder builder; 
  builder.BuildWorld();
  builder.BuildHead();
  builder.BuildTrunk();
  return builder.GetPhantom(); 
}


void  G4HumanPhantomConstruction::SetBodyPartSensitivity(G4String newBodyPartLV, G4bool bodypartSensitivity)
{

}
