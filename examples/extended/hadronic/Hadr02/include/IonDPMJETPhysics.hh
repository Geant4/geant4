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
// $Id: IonDPMJETPhysics.hh,v 1.0 2010/08/26 10:51:25 antoni Exp $
// GRAS tag $Name: gras-02-05-02 $
//
//---------------------------------------------------------------------------
//
// Header:    IonDPMJETPhysics
//
// Author:    copy from P.Truscott manuel DPMJET2.5 
//
// 
// Customer:          
// Contract:          
//
// Modifications are provided according to
//
// Organisation:        
// Customer:            
// Contract:            
//
// Modified:     26.08.2010
//
// ------------------------------------------------------------
//

#ifndef IonDPMJETPhysics_h
#define IonDPMJETPhysics_h 1

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4BinaryLightIonReaction;
class G4DPMJET2_5Model;
class G4VCrossSectionDataSet;

class IonDPMJETPhysics : public G4VHadronPhysics
{
public:
  IonDPMJETPhysics(G4bool val);
  virtual ~IonDPMJETPhysics();

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
  G4DPMJET2_5Model*          theDPM;
  G4bool isBinary;
};

#endif
