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

#ifndef G4VFieldPropagation_h
#define  G4VFieldPropagation_h 1

#include "G4KineticTrackVector.hh"
#include "G4V3DNucleus.hh"

class G4VFieldPropagation
{
public:
  G4VFieldPropagation();
  virtual ~G4VFieldPropagation();

private:
  G4VFieldPropagation(const  G4VFieldPropagation &right);
  const G4VFieldPropagation & operator=(const G4VFieldPropagation & right);
  G4int operator==(const G4VFieldPropagation & right) const;
  G4int operator!=(const G4VFieldPropagation & right) const;

public:
  virtual void Init(G4V3DNucleus * theNucleus) = 0;
  virtual void Transport(G4KineticTrackVector &theActive,
			 const G4KineticTrackVector &theSpectators,
			 G4double theTimeStep) = 0;
  virtual G4ThreeVector GetMomentumTransfer() const =0;
};

#endif

