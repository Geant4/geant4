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
// $Id: G4FissionBarrier.hh 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up
// 21.03.2013 V.Ivanchenko redesign parameters classes

#ifndef G4FissionBarrier_h
#define G4FissionBarrier_h 1

#include "G4VFissionBarrier.hh"
#include "globals.hh"
#include "G4CameronShellPlusPairingCorrections.hh"

class G4FissionBarrier : public G4VFissionBarrier
{
public:

  explicit G4FissionBarrier();

  virtual ~G4FissionBarrier();

  virtual G4double FissionBarrier(G4int A, G4int Z, G4double U) const final;

private:

  G4double BarashenkovFissionBarrier(G4int A, G4int Z) const;
  
  G4FissionBarrier(const G4FissionBarrier & right) = delete;
  const G4FissionBarrier & operator=(const G4FissionBarrier & right) = delete;
  G4bool operator==(const G4FissionBarrier & right) const = delete;
  G4bool operator!=(const G4FissionBarrier & right) const = delete;

  const G4CameronShellPlusPairingCorrections* SPtr;
};

#endif
