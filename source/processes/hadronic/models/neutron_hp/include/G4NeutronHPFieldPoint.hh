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
// $Id: G4NeutronHPFieldPoint.hh,v 1.5 2001-07-26 09:27:59 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFieldPoint_h
#define G4NeutronHPFieldPoint_h 1

#include "globals.hh"

class G4NeutronHPFieldPoint
{
  public:
  
  G4NeutronHPFieldPoint()
  {
    X = 0;
    nP = 0;
    Y = NULL;
  }
  
  G4NeutronHPFieldPoint(G4int n);
  
  void operator= (const G4NeutronHPFieldPoint & aSet);

  ~G4NeutronHPFieldPoint();
    
  void InitY(G4int n);

  inline G4int GetDepth() const {return nP;}
  inline G4double GetX() const {return X;}
  inline G4double GetY(G4int i) const {return Y[i];}
  
  inline void SetX(G4double e) {X = e;}
  inline void SetY(G4int i, G4double x) {Y[i] = x;}
  
  inline void SetData(G4double e, G4int i, G4double x) {X = e; Y[i] = x;}
  
  private:
  
  G4double X;
  G4double * Y;
  G4int nP;
};

#endif
