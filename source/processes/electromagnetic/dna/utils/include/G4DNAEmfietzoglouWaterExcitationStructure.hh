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
// Based on the work described in
// Rad Res 163, 98-111 (2005)
// D. Emfietzoglou, H. Nikjoo
// 
// Authors of the class (2014):
// I. Kyriakou (kyriak@cc.uoi.gr)
// D. Emfietzoglou (demfietz@cc.uoi.gr)
// S. Incerti (incerti@cenbg.in2p3.fr)
//

#ifndef G4DNAEmfietzoglouWaterExcitationStructure_h
#define G4DNAEmfietzoglouWaterExcitationStructure_h 1
 
#include "globals.hh"
#include <vector>

class G4DNAEmfietzoglouWaterExcitationStructure
{

public:
  
  G4DNAEmfietzoglouWaterExcitationStructure();
  
  virtual ~G4DNAEmfietzoglouWaterExcitationStructure();
  
  G4double ExcitationEnergy(G4int level);

  G4int NumberOfLevels() { return nLevels; }
      
private:
   
  G4int nLevels;

  std::vector<G4double> energyConstant;
  
};

#endif
