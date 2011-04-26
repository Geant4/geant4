
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software.                                *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This code implementation is the intellectual property of the ESA.*
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4IonFTFPBinaryCascadePhysics.hh,v 1.0 2010/08/26 10:51:25 antoni Exp $
// GRAS tag $Name: $
//
//---------------------------------------------------------------------------
//
// Header:    G4IonFTFPBinaryCascadePhysics
//
// Author:    V.Ivanchenko  02.03.2011
//
// 
//
// Modified:     
//
// ------------------------------------------------------------
//

#ifndef G4IonFTFPBinaryCascadePhysics_h
#define G4IonFTFPBinaryCascadePhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4FTFBuilder;
class G4BinaryLightIonReaction;

class G4IonFTFPBinaryCascadePhysics : public G4VPhysicsConstructor
{
public:

  G4IonFTFPBinaryCascadePhysics(G4int ver = 0);
  virtual ~G4IonFTFPBinaryCascadePhysics();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  void ConstructParticle();
  void ConstructProcess();

private:

  void AddProcess(const G4String&, G4ParticleDefinition*, G4bool isIon);

  G4VCrossSectionDataSet*   fTripathi;
  G4VCrossSectionDataSet*   fTripathiLight;
  G4VCrossSectionDataSet*   fShen;
  G4VCrossSectionDataSet*   fIonH;
  G4BinaryLightIonReaction* theIonBC;
  G4HadronicInteraction*    theFTFP;
  G4FTFBuilder*             theBuilder;

  G4int  verbose;
  G4bool wasActivated;
};


#endif








