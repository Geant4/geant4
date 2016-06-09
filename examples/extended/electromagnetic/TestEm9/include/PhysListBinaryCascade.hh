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
// $Id: PhysListBinaryCascade.hh,v 1.5 2003/12/05 11:22:45 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ------------------------------------------------------------
//	GEANT 4 class header file
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//      It is provide PhysicsList for Binary Cascade for
//      protons and neutrons with the energy E<3 GeV
//
// ------------------------------------------------------------
//	History
//        Created:       14.10.02  V.Ivanchenko
//
//        Modified:
//
// ------------------------------------------------------------
//
#ifndef PhysListBinaryCascade_h
#define PhysListBinaryCascade_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

class PhysListBinaryCascade : public G4VPhysicsConstructor
{
  public:
    PhysListBinaryCascade(const G4String& name = "binary");
    virtual ~PhysListBinaryCascade();

  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    void ConstructParticle() {};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess();

  private:

    G4ProtonInelasticProcess       theIPproton;
    G4NeutronInelasticProcess      theIPneutron;
    G4ProtonInelasticCrossSection  thePXSec;
    G4NeutronInelasticCrossSection theNXSec;

};

#endif
