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

#ifndef G4VPrimitivePlotter_H
#define G4VPrimitivePlotter_H 1

#include "G4VPrimitiveScorer.hh"
#include <map>

class G4VPrimitivePlotter : public G4VPrimitiveScorer
{
 public:  // with description
  G4VPrimitivePlotter(G4String name, G4int depth = 0)
    : G4VPrimitiveScorer(name, depth)
  {
    ;
  }
  virtual ~G4VPrimitivePlotter() { ; }

 public:
  void Plot(G4int copyNo, G4int histID) { hitIDMap[copyNo] = histID; }

 protected:
  std::map<G4int, G4int> hitIDMap;

 public:
  G4int GetNumberOfHist() const { return (G4int)hitIDMap.size(); }
};

#endif
