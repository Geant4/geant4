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
// $Id: G4WeightCutOffProcess.hh,v 1.2 2002-10-16 16:26:59 dressel Exp $
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

class G4VGCellFinder;
class G4VIStore;

class G4WeightCutOffProcess : public G4VProcess
{

public:  // with description

  G4WeightCutOffProcess(G4double wsurvival,
			G4double wlimit,
			G4double isource,
			G4VIStore *istore,
			const G4VGCellFinder &aGCellFinder,
			const G4String &aName = "WeightCutOffProcess");
    // create a G4ParticleChange

  virtual ~G4WeightCutOffProcess();
    // delete the G4ParticleChange

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
    // make the process beeing forced

  virtual G4VParticleChange * PostStepDoIt(const G4Track&, 
				   const G4Step&);


  const G4String &GetName() const;

 
public:  // without description

  //  no operation in  AtRestDoIt and  AlongStepDoIt

  virtual G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*);
  
  virtual G4double 
  AtRestGetPhysicalInteractionLength(const G4Track& ,
				     G4ForceCondition*);

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);
  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);
  
private:

  G4WeightCutOffProcess(const G4WeightCutOffProcess &);
  G4WeightCutOffProcess &operator=(const G4WeightCutOffProcess &);

  G4ParticleChange *fParticleChange;
  G4double fWeightSurvival;
  G4double fWeightLimit;
  G4double fSourceImportance;
  G4VIStore *fIStore;
  const G4VGCellFinder &fGCellFinder;

};

#endif





