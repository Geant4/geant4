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
#ifndef G4BertiniElasticCollision_h
#define G4BertiniElasticCollision_h 1

#include "G4BertiniModel.hh"

class G4BertiniElasticCollision : public G4BertiniModel {

public:
  G4BertiniElasticCollision();
  ~G4BertiniElasticCollision();
  void interpolateElasticNeutronData(G4int medium, G4int kdd, G4double e);
  void geti(G4double es, G4int npts, G4double e, G4int i);
  void scatteringWithHydrogen(); 

private:
  G4double pt[16];
  G4double col[23];
};

#endif

  












