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
// $Id: PhysListGN.hh,v 1.1 2004-01-07 11:30:00 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//	GEANT 4 class header file
//
//      This class is an derived class of G4VPhysicsConstructor
//      It is provide PhysicsList for gamma nuclear interactions
//      for the energy E<3.6 GeV
//
//        Created:       05.12.03  V.Ivanchenko
//
//        Modified:
//
// ------------------------------------------------------------
//
#ifndef PhysListGN_h
#define PhysListGN_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4PhotoNuclearProcess.hh"

class G4GammaNuclearReaction;

class PhysListGN : public G4VPhysicsConstructor
{
  public:
    PhysListGN(const G4String& name = "photoNuc");
    virtual ~PhysListGN();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    void ConstructParticle() {};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess();

  private:

    G4PhotoNuclearProcess thePhotoNuclearProcess;
    G4GammaNuclearReaction * theGammaReaction;
};

#endif
