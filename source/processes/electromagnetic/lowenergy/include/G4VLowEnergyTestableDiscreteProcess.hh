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
// $Id: G4VLowEnergyTestableDiscreteProcess.hh,v 1.2 2006-05-25 17:57:10 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VLowEnergyTestableDiscreteProcess_hh
#define G4VLowEnergyTestableDiscreteProcess_hh 1
 
#include "G4VDiscreteProcess.hh"
 
class G4VLowEnergyTestableDiscreteProcess : public G4VDiscreteProcess
{
public:

  G4VLowEnergyTestableDiscreteProcess(const G4String & name) : G4VDiscreteProcess(name) {}

  virtual ~G4VLowEnergyTestableDiscreteProcess() {}
 
  G4double DumpMeanFreePath(const G4Track & aTrack, G4double previousStepSize, G4ForceCondition * condition)
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }
   
private:
  // Hides default constructor and assignment operator as private 
  G4VLowEnergyTestableDiscreteProcess();
  G4VLowEnergyTestableDiscreteProcess & operator=(const G4VLowEnergyTestableDiscreteProcess & right);
};

#endif /* G4VLowEnergyTestableDiscreteProcess_hh */

