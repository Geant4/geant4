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
//    **************************************
//    *                                    *
//    *    RemSimHadronicBinary.hh        *
//    *                                    *
//    **************************************
//
// $Id: RemSimHadronicBinary.hh,v 1.5 2006-11-15 18:39:30 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 

#ifndef RemSimHadronicBinary_h
#define RemSimHadronicBinary_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4ProtonInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4QGSModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

class RemSimHadronicBinary: public G4VPhysicsConstructor 
{
  public:
    RemSimHadronicBinary(const G4String& name = "hadronic-binary");
    virtual ~RemSimHadronicBinary();

  protected:
    // Construct particle and physics
    void ConstructParticle(){};
    void ConstructProcess();

 private:
   G4ProtonInelasticCrossSection protonCrossSection;
   G4PiNuclearCrossSection pionCrossSection;
   G4NeutronInelasticCrossSection neutronCrossSection;
   G4TheoFSGenerator* QGSPModel;
   G4GeneratorPrecompoundInterface* theCascade;
   G4ExcitationHandler theHandler;
   G4PreCompoundModel* thePreEquilib;
   G4QGSMFragmentation* theFragmentation;
   G4ExcitedStringDecay* theStringDecay;
   G4QGSModel<G4QGSParticipants> theStringModel;
};
#endif








