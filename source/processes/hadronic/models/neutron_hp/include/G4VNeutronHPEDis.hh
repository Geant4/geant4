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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VNeutronHPEDis.hh,v 1.6 2001-07-26 09:28:21 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4VNeutronHPEDis_h
#define G4VNeutronHPEDis_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"

class G4VNeutronHPEDis
{
  public:
  G4VNeutronHPEDis()
  {
  }
  virtual ~G4VNeutronHPEDis()
  {
  }
  
  virtual void Init(G4std::ifstream & theData) = 0;
  
  virtual G4double GetFractionalProbability(G4double anEnergy) = 0;
  
  virtual G4double Sample(G4double anEnergy) = 0; 
  
  private:
    
};

#endif
