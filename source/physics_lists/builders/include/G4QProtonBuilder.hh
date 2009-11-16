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
// ClassName:   G4QProtonBuilder
//
// Author: 2009 M. Kosov
//
// Modified:
//
//-----------------------------------------------------------------------------
// Short description: for G4QDiscProcessMixer use in the QGSC_CHIPS physics list
//-----------------------------------------------------------------------------
//
#ifndef G4QProtonBuilder_h
#define G4QProtonBuilder_h 1

#include "globals.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4QInelastic.hh"
#include "G4QDiscProcessMixer.hh"
#include "G4VProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include <vector>

class G4QProtonBuilder
{
public:
  G4QProtonBuilder();
  virtual ~G4QProtonBuilder();

public: 
  void Build();
  void RegisterMe(G4VProtonBuilder * aB) {theModelCollections.push_back(aB);}

private:
  G4QDiscProcessMixer*            theProcessMixer;
  std::vector<G4VProtonBuilder *> theModelCollections;
  G4ProtonInelasticProcess*       theProtonInelastic;
  G4QInelastic*                   theCHIPSInelastic;    

  G4bool wasActivated;
};

// 2009 by M. Kosov

#endif
