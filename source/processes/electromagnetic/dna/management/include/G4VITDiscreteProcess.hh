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
/*
 * G4VITDiscreteProcess.hh
 *
 *  Created on: 9 août 2015
 *      Author: matkara
 */

#ifndef SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITDISCRETEPROCESS_HH_
#define SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITDISCRETEPROCESS_HH_

#include <G4VITProcess.hh>
#include "globals.hh"
#include "G4ios.hh"

class G4VITDiscreteProcess: public G4VITProcess
{
  //  Abstract class which defines the public behavior of
  //  discrete physics interactions.
public:

  G4VITDiscreteProcess(const G4String&, G4ProcessType aType = fNotDefined);
  G4VITDiscreteProcess(G4VITDiscreteProcess &);

  virtual ~G4VITDiscreteProcess();

public:
  // with description
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                        G4double previousStepSize,
                                                        G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                         G4double,
                                                         G4double,
                                                         G4double&,
                                                         G4GPILSelection*)
  {
    return -1.0;
  }
  ;

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&,
                                                      G4ForceCondition*)
  {
    return -1.0;
  }
  ;

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&)
  {
    return 0;
  }
  ;

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&)
  {
    return 0;
  }
  ;

protected:
  // with description
  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition)=0;
  //  Calculates from the macroscopic cross section a mean
  //  free path, the value is returned in units of distance.

private:
  // hide default constructor and assignment operator as private
  G4VITDiscreteProcess();
  G4VITDiscreteProcess & operator=(const G4VITDiscreteProcess &right);

};

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITDISCRETEPROCESS_HH_ */
