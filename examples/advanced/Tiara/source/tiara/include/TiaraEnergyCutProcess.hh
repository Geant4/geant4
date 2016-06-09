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
// $Id: TiaraEnergyCutProcess.hh,v 1.4 2006/06/29 15:43:37 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
