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
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//
// ------------------------------------------------------------ 
//	History
//        Created:       14.10.02  V.Ivanchenko
//
//        Modified:
// 
// ------------------------------------------------------------
//
#ifndef test31GeneralPhysics_h
#define test31GeneralPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"
#include "globals.hh"

class test31GeneralPhysics : public G4VPhysicsConstructor
{
  public: 
    test31GeneralPhysics(const G4String& name = "general", 
                         G4VProcess* test31StepCut = 0);
    virtual ~test31GeneralPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  private:
    G4Decay fDecayProcess;
    G4VProcess* stepCut;
};


#endif








