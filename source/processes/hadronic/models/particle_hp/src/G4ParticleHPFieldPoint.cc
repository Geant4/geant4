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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//

#include "G4ParticleHPFieldPoint.hh"

G4ParticleHPFieldPoint::G4ParticleHPFieldPoint(G4int n)
  {
    nP = n;
    X = 0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
  
void G4ParticleHPFieldPoint::operator= (const G4ParticleHPFieldPoint & aSet)
  {
    if(&aSet!=this)
    {
      X = aSet.GetX();
      delete [] Y;
      Y = new G4double[aSet.GetDepth()];
      for(G4int i=0; i<aSet.GetDepth(); i++) Y[i] = aSet.GetY(i);
    }
  }

G4ParticleHPFieldPoint::~G4ParticleHPFieldPoint()
  {
    delete [] Y;
  }
    
void G4ParticleHPFieldPoint::InitY(G4int n)
  {
    nP = n;
    X=0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
