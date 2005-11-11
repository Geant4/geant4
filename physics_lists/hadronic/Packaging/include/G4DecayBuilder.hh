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
// $Id: G4DecayBuilder.hh,v 1.1 2005-11-11 22:56:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4DecayBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
//
//----------------------------------------------------------------------------
//

#ifndef G4DecayBuilder_h
#define G4DecayBuilder_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class G4Decay;

class G4DecayBuilder : public G4VPhysicsConstructor
{
  public: 
    G4DecayBuilder(const G4String& name = "decay");
    virtual ~G4DecayBuilder();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
  virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:
  G4Decay* fDecayProcess;
  G4bool   wasActivated;
};

#endif








