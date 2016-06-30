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
// $Id: G4VFissionBarrier.hh 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#ifndef G4VFissionBarrier_h
#define G4VFissionBarrier_h 1

#include "globals.hh"
#include "G4HadronicException.hh"


class G4VFissionBarrier
{
public:
  G4VFissionBarrier();
  virtual ~G4VFissionBarrier();

private:
  G4VFissionBarrier(const G4VFissionBarrier & right) = delete;

  const G4VFissionBarrier & operator=(const G4VFissionBarrier & right) = delete;
  G4bool operator==(const G4VFissionBarrier & right) const = delete;
  G4bool operator!=(const G4VFissionBarrier & right) const = delete;
  
public:

  virtual G4double FissionBarrier(G4int A, G4int Z, G4double U) const = 0;
  

};

#endif
