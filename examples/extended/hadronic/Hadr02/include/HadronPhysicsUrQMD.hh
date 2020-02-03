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
/// \file hadronic/Hadr02/include/HadronPhysicsUrQMD.hh
/// \brief Definition of the HadronPhysicsUrQMD class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2012 Andrea Dotti
//   created from HadronPhysicsFTFP_BERT
// Modified:
// 07.02.2012 A. Dotti: First version
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsUrQMD_h
#define HadronPhysicsUrQMD_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4PiKBuilder.hh"
#include "UrQMDPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "UrQMDProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "UrQMDNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"

#include "G4AntiBarionBuilder.hh"
#include "UrQMDAntiBarionBuilder.hh"

class HadronPhysicsUrQMD : public G4VPhysicsConstructor
{
public: 

  HadronPhysicsUrQMD(G4int verbose =1);
  virtual ~HadronPhysicsUrQMD();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:

  void CreateModels();
  G4HadronicProcess* FindInelasticProcess(const G4ParticleDefinition*);
    
  G4NeutronBuilder * fNeutrons;
  UrQMDNeutronBuilder * fUrQMDNeutron;
 
  G4PiKBuilder * fPiK;
  UrQMDPiKBuilder * fUrQMDPiK;
    
  G4ProtonBuilder * fPro;
  UrQMDProtonBuilder * fUrQMDPro;    
    
  G4HyperonFTFPBuilder * fHyperon;
    
  G4AntiBarionBuilder * fAntiBaryon;
  UrQMDAntiBarionBuilder * fUrQMDAntiBaryon;

};

#endif

