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
// $Id: Test2RunSD.hh,v 1.1 2010-11-03 08:48:57 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Test2RunSD_h
#define Test2RunSD_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "Test2PhantomHit.hh"
#include "Test2SDHitSum.hh"

class G4Event;

class Test2RunSD : public G4Run
{

public:
  Test2RunSD(const G4String& detName,const G4String& hcname,
	     std::vector<G4String>& hcnameVec);
  virtual ~Test2RunSD();

public:
  virtual void RecordEvent(const G4Event*);
  void DumpQuantitiesToFile();
  G4double GetTotal(G4int i) const;

private:
  //-- traditional sensitive detetector
  G4int fSdID;
  //
  Test2SDHitSum* HitSum;

};

#endif

