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
//
//    ************************
//    *                      *
//    *    HadrontherapyDecay.hh    *
//    *                      *           
//    ************************
//
// $Id: HadrontherapyDecay.hh,v 1.3 2005-05-18 07:53:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifndef HADRONTHERAPYDECAY_HH
#define HADRONTHERAPYDECAY_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"
#include "G4Decay.hh"

class HadrontherapyDecay : public G4VPhysicsConstructor {

public: 

  HadrontherapyDecay(const G4String& name = "decay");
  
  virtual ~HadrontherapyDecay();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();  
};
#endif

