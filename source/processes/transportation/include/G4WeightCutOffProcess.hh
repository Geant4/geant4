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
// $Id: G4WeightCutOffProcess.hh,v 1.1 2002-10-10 13:25:31 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffProcess
//
// Class description:
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightCutOffProcess_hh 
#define G4WeightCutOffProcess_hh G4WeightCutOffProcess_hh

#include "G4VProcess.hh"
#include "G4VTrackTerminator.hh"
#include "G4GeometryCell.hh"

class G4VParallelStepper;
class G4VIStore;

class G4WeightCutOffProcess : public G4VProcess
{

public:  // with description

  G4WeightCutOffProcess(G4double wsurvival,
			G4double wlimit,
			G4double isource,
			G4VIStore *istore,
			G4VParallelStepper  *astepper,
			const G4String &aName = "WeightCutOffProcess");
    // create a G4ParticleChange

  ~G4WeightCutOffProcess();
    // delete the G4ParticleChange

  G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
    // make the process beeing forced

  G4VParticleChange * PostStepDoIt(const G4Track&, 
				   const G4Step&);


  G4String GetName() const {
    return theProcessName;
  }

 
public:  // without description

  //  no operation in  AtRestDoIt and  AlongStepDoIt

  G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                             G4double  ,
                                             G4double  ,
                                             G4double& ,
                                             G4GPILSelection*) {return -1.0;}
  
  G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
                                              G4ForceCondition*) {return -1.0;}
  
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) {return 0;}
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) {return 0;}
  
private:

  G4WeightCutOffProcess(const G4WeightCutOffProcess &);
  G4WeightCutOffProcess &operator=(const G4WeightCutOffProcess &);

  G4GeometryCell GetPostGeometryCell(const G4Step &aStep);

  G4double fWeightSurvival;
  G4double fWeightLimit;
  G4double fSourceImportance;
  G4VIStore *fIStore;
  G4VParallelStepper *fPstepper;
  G4ParticleChange *fParticleChange;

};

#endif





