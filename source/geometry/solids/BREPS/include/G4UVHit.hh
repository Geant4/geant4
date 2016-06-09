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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4UVHit
//
// Class Description:
//   
// Utility class for knots.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4UV_Hit
#define __G4UV_Hit

#include "globals.hh"

class G4UVHit
{
public:
    G4UVHit() {u=-1; v=-1; next=0;}	
    G4UVHit(G4double u_hit, G4double v_hit) {u = u_hit; v = v_hit; next=0;}
    ~G4UVHit() {}

    void SetNext(G4UVHit* n) { next = n; }
    G4UVHit* GetNext() { return next; }
    G4double GetU() const { return u; }
    G4double GetV() const { return v; }
    
private:

    G4UVHit* next;
    G4double u, v;
};

#endif
