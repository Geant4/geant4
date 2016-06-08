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

  












