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
// $Id: ExN07SteppingVerbose.hh,v 1.1 2006-11-04 19:23:07 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

class ExN07SteppingVerbose;

#ifndef ExN07SteppingVerbose_h
#define ExN07SteppingVerbose_h 1

#include <vector>
#include "G4VSteppingVerbose.hh"
#include "G4SliceTimer.hh"

class G4Region;

class ExN07SteppingVerbose : public G4VSteppingVerbose
{
 public:   

   ExN07SteppingVerbose();
  ~ExN07SteppingVerbose();

  void InitializeTimers();
  void Report();

  virtual void NewStep();
  virtual void StepInfo();

  // Following methods are not used
  virtual void TrackBanner();
  virtual void AtRestDoItInvoked();
  virtual void AlongStepDoItAllDone();
  virtual void PostStepDoItAllDone();
  virtual void AlongStepDoItOneByOne();
  virtual void PostStepDoItOneByOne();
  virtual void TrackingStarted();
  virtual void DPSLStarted();
  virtual void DPSLUserLimit();
  virtual void DPSLPostStep();
  virtual void DPSLAlongStep();
  virtual void VerboseTrack();
  virtual void VerboseParticleChange();

 private:
  G4int FindRegion(G4Region*);

 private:
  std::vector<G4SliceTimer*> fTimers;
  G4int nRegions,nTimers,regIdx;
  G4bool ep;
};


#endif
