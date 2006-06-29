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
//
// $Id: G4ProcessTest.hh,v 1.6 2006-06-29 19:48:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test of electromagnetic physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4PROCESSTEST_HH
#define G4PROCESSTEST_HH 1

#include "globals.hh"

class G4VProcess;
class G4ProcessManager;
class G4ParticleDefinition;
class G4Track;
class G4Step;

class G4ProcessTest {
 
public:

  G4ProcessTest();

  virtual ~G4ProcessTest();
 
  void buildTables(const G4String& type, G4bool isPolarised = false);

  void postStepTest(const G4Track& track,const G4Step& step) const;

  void alongStepTest(const G4Track& track, const G4Step& step) const;

  G4ParticleDefinition* getIncidentParticleType() { return def; }

  protected:

  virtual G4VProcess* createProcess() = 0;
  virtual G4VProcess* createBremsstrahlung() = 0;
  virtual G4VProcess* createElectronIonisation() = 0;
  virtual G4ParticleDefinition* createIncidentParticle();

private:
  
  // Hide copy constructor and assignment operator
  G4ProcessTest(const G4ProcessTest&);
  G4ProcessTest & operator=(const G4ProcessTest &right);

  G4VProcess* process;
  G4VProcess* ioni;
  G4VProcess* brem;
  G4ProcessManager* eProcessManager;
  G4ProcessManager* gProcessManager;
  G4ParticleDefinition* def;      // not owned pointer
};
 
#endif

