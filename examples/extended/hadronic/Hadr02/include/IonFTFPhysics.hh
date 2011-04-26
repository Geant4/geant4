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
// $Id: IonFTFPhysics.hh,v 1.0 2010/08/26 10:51:25 antoni Exp $
// GRAS tag $Name: gras-02-05-02 $
//
//---------------------------------------------------------------------------
//
// Header:    IonFTFPhysics
//
// Author:    V.Ivanchenko  02.03.2011
//
// 
//
// Modified:     
//
// ------------------------------------------------------------
//

#ifndef IonFTFPhysics_h
#define IonFTFPhysics_h 1

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4VCrossSectionDataSet;
class G4BinaryLightIonReaction;

class IonFTFPhysics : public G4VHadronPhysics
{
public:
  IonFTFPhysics(G4bool val);
  virtual ~IonFTFPhysics();

public:
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  void ConstructProcess();

private:

  G4VCrossSectionDataSet* fTripathi;
  G4VCrossSectionDataSet* fTripathiLight;
  G4VCrossSectionDataSet* fShen;
  G4VCrossSectionDataSet* fIonH;
  G4BinaryLightIonReaction*  theIonBC;
  G4bool isBinary;
};


#endif








