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
// $Id: G4ScoreLogColorMap.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4ScoreLogColorMap_h
#define G4ScoreLogColorMap_h 1

#include "globals.hh"
#include "G4VScoreColorMap.hh"

class G4ScoreLogColorMap : public G4VScoreColorMap
{
  public:
      G4ScoreLogColorMap(G4String mName);
      virtual ~G4ScoreLogColorMap();

  public:
      virtual void GetMapColor(G4double val, G4double color[4]);

  // draw a color chart
  virtual void DrawColorChartBar(G4int nPoint);
  virtual void DrawColorChartText(G4int nPoint);

};

#endif

