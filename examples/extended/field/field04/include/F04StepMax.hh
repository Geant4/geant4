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
// $Id: F04StepMax.hh 68021 2013-03-13 13:36:07Z gcosmo $
//
/// \file field/field04/include/F04StepMax.hh
/// \brief Definition of the F04StepMax class
//

#ifndef F04StepMax_h
#define F04StepMax_h 1

#include "globals.hh"

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"

class F04StepMax : public G4VDiscreteProcess
{
  public:

    F04StepMax(const G4String& processName = "UserStepMax");
    F04StepMax(F04StepMax &);

    virtual ~F04StepMax();

    virtual G4bool IsApplicable(const G4ParticleDefinition&);

    void SetStepMax(G4double);

    G4double GetStepMax() {return fMaxChargedStep;};

    virtual G4double PostStepGetPhysicalInteractionLength(const G4Track&,
                                                          G4double,
                                                          G4ForceCondition*);

    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  protected:

    virtual G4double GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition*);

  private:

    // hide assignment operator as private
    F04StepMax & operator=(const F04StepMax &right);
    F04StepMax(const F04StepMax&);

  private:

    G4double fMaxChargedStep;

};

#endif
