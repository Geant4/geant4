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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QNeutronHPBuilder
//
// Author: April 2012 M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
// Short description: for use in CHIPS_HP physics list (mix Q w/ HP at 20 MeV)
//----------------------------------------------------------------------------
//
#ifndef G4QNeutronHPBuilder_h
#define G4QNeutronHPBuilder_h 1

#include "globals.hh"

#include "G4QInelastic.hh"
#include "G4QNGamma.hh"
//#include "G4QFission.hh"
#include "G4QDiscProcessMixer.hh"
#include "G4NeutronBuilder.hh"
#include "G4VNeutronBuilder.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include <vector>

class G4QNeutronHPBuilder
{
public: 
  G4QNeutronHPBuilder();
  virtual ~G4QNeutronHPBuilder();

public: 
  void Build();
  void RegisterMe(G4VNeutronBuilder * aB) {theModelCollections.push_back(aB);}

private:
  G4QDiscProcessMixer*            theInProcessMixer;
  G4QDiscProcessMixer*            theNgProcessMixer;
  G4QDiscProcessMixer*            theFiProcessMixer;
  G4NeutronBuilder*               theNeutrons;
  std::vector<G4VNeutronBuilder*> theModelCollections;
  G4NeutronInelasticProcess*      theNeutronInelastic;
  G4HadronFissionProcess*         theNeutronFission;
  G4HadronCaptureProcess*         theNeutronCapture;
  G4QInelastic*                   theCHIPSInelastic;
  G4QNGamma*                      theCHIPSNGamma;
  //G4QFission*                     theCHIPSFission;
  G4NeutronHPBuilder*             theHPNeutron;
  
  G4bool wasActivated;

};

// 2012 by M. Kosov

#endif

