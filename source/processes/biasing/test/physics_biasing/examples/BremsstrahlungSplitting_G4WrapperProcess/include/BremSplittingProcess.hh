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
// $Id: BremSplittingProcess.hh,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, May 2007
//
#ifndef BREMSPLITTINGPROCESS_HH
#define BREMSPLITTINGPROCESS_HH 

#include "G4WrapperProcess.hh"

struct BremSplittingProcess : public G4WrapperProcess {

  BremSplittingProcess();	
  
  virtual ~BremSplittingProcess();	
  
  // Override PostStepDoIt  method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

};

#endif
