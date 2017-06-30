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
// $Id: G4EmParticleList 99930 2016-10-11 16:34:52Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmParticleList 
//
// Author:      V.Ivanchenko 04.05.2017 
//
// Modified:
//
//----------------------------------------------------------------------------
//
// class description
//
//  The class provides list of particles which have EM interactions
//
// class description - end


#ifndef G4EmParticleList_h
#define G4EmParticleList_h 1

#include "globals.hh"
#include <vector>

class G4EmParticleList {

public:

  explicit G4EmParticleList();

  ~G4EmParticleList();

  const std::vector<G4String>& PartNames() const;

private:
  std::vector<G4String>  pNames; 
  
};

#endif

