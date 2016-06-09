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
// $Id: TiaraEnergyCutProcess.hh,v 1.3 2003/06/25 09:12:39 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraEnergyCutProcess
//

#ifndef TiaraEnergyCutProcess_hh 
#define TiaraEnergyCutProcess_hh TiaraEnergyCutProcess_hh 

#include "G4VProcess.hh"
#include "G4VTrackTerminator.hh"

class G4VScorer;

class TiaraEnergyCutProcess : public G4VProcess
{

public:  // with description

  explicit TiaraEnergyCutProcess(G4double minEnergyCut);

  virtual ~TiaraEnergyCutProcess();

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
    // make processed being forced

  virtual G4VParticleChange * PostStepDoIt(const G4Track&, const G4Step&);
    // message "scorer" with  G4Step and a G4GeometryCellStep from the "mass" 
    // geometry


    // to be called by the importance process if the track should
    // be killed after scoring
  virtual const G4String &GetName() const ;

public:  // without description

  // no operation in  AtRestDoIt and  AlongStepDoIt

  virtual G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*);
  
  virtual G4double 
  AtRestGetPhysicalInteractionLength(const G4Track&,
				     G4ForceCondition*);
  
  virtual G4VParticleChange* AtRestDoIt(const G4Track&,
					const G4Step&);
  virtual G4VParticleChange* AlongStepDoIt(const G4Track&,
					   const G4Step&);
  
private:

  TiaraEnergyCutProcess(const TiaraEnergyCutProcess &);
  TiaraEnergyCutProcess &operator=(const TiaraEnergyCutProcess &);
  
private:

  G4double fMinEnergyCut;
};

#endif
