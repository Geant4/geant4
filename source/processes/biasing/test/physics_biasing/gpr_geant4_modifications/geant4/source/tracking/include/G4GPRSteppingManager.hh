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
// $Id: G4GPRSteppingManager.hh,v 1.1 2007-09-06 22:19:22 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//---------------------------------------------------------------
//
// G4SteppingManager.hh
//
// class description:
//  This is the class which plays an essential role in tracking 
//  the particle. It takes cares all message passing between
//  objects in the different categories (for example, 
//  geometry(including transportation), interactions in 
//  matter, etc). It's public method 'stepping' steers to step 
//  the particle.
//  Geant4 kernel use only
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige


class G4GPRSteppingManager;

#ifndef G4GPRSteppingManager_h
#define G4GPRSteppingManager_h 1

#include "G4SteppingManager.hh"
#include <vector>
#include "G4GPRProcessWrappers.hh"

class G4GPRSteppingManager : public G4SteppingManager {

public: 

  G4GPRSteppingManager() {};
  
  ~G4GPRSteppingManager() {};

   virtual G4StepStatus Stepping();

protected:

  virtual void DefinePhysicalStepLength();
  virtual void InvokeAtRestDoItProcs();
  virtual void InvokeAlongStepDoItProcs();
  virtual void InvokePostStepDoItProcs();
  virtual void InvokePSDIP(size_t np);

private:
  std::vector<G4GPRProcessWrappers::G4GPRDiscreteGPIL>* pDiscreteGPIL;
  std::vector<G4GPRProcessWrappers::G4GPRContinuousGPIL>* pContinuousGPIL;
  std::vector<G4GPRProcessWrappers::G4GPRAtRestGPIL>* pAtRestGPIL;

  std::vector<G4GPRProcessWrappers::G4GPRAtRestDoIt>* pAtRestDoIt;
  std::vector<G4GPRProcessWrappers::G4GPRContinuousDoIt>* pContinuousDoIt;
  std::vector<G4GPRProcessWrappers::G4GPRDiscreteDoIt>* pDiscreteDoIt;
};

#endif
