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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//

#include "G4NeutronHPFieldPoint.hh"

G4NeutronHPFieldPoint::G4NeutronHPFieldPoint(G4int n)
  {
    nP = n;
    X = 0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
  
void G4NeutronHPFieldPoint::operator= (const G4NeutronHPFieldPoint & aSet)
  {
    if(&aSet!=this)
    {
      X = aSet.GetX();
      if(Y!=NULL) delete [] Y;
      Y = new G4double[aSet.GetDepth()];
      for(G4int i=0; i<aSet.GetDepth(); i++) Y[i] = aSet.GetY(i);
    }
  }

G4NeutronHPFieldPoint::~G4NeutronHPFieldPoint()
  {
   if(Y!=NULL) delete [] Y;
  }
    
void G4NeutronHPFieldPoint::InitY(G4int n)
  {
    nP = n;
    X=0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
