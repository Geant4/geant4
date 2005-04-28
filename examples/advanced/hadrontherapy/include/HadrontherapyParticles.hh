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
// $Id: HadrontherapyParticles.hh,v 1.2 2005-04-28 20:39:33 mpiergen Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// -------------------------------------------------------------------

#ifndef HadrontherapyPARTICLES_HH
#define HadrontherapyPARTICLES_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class HadrontherapyParticles : public G4VPhysicsConstructor {
public: 
  HadrontherapyParticles(const G4String& name = "particles");
  
  virtual ~HadrontherapyParticles();
  
  virtual void ConstructParticle();
  
  // This method is dummy
  virtual void ConstructProcess() {};
  
};

#endif








