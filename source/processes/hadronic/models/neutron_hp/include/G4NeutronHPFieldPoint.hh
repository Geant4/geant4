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
// $Id: G4NeutronHPFieldPoint.hh,v 1.8 2006-06-29 20:47:43 gunter Exp $
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
