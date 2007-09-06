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
// $Id: G4GPRMask.hh,v 1.2 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRMASK_HH
#define G4GPRMASK_HH

class G4GPRMask {

public:

  G4GPRMask(const G4String& name, std::vector<unsigned> processIndices)
    :fName(name)
    ,fActive(true) 
    ,fProcessIndices(processIndices)
  {}

  G4bool IsActive() {return fActive;}

  void ChangeState() 
  {
    fActive = !fActive;
  }

  const std::vector<unsigned>& GetProcessIndices() 
  {
    return fProcessIndices;
  }

  G4String GetName() {return fName;}

private:

  G4String fName;
  G4bool fActive;
  std::vector<unsigned> fProcessIndices;
};

#endif
