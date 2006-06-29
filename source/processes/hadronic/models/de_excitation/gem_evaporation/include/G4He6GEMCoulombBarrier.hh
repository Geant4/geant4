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
// $Id: G4He6GEMCoulombBarrier.hh,v 1.3 2006-06-29 20:16:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#ifndef G4He6GEMCoulombBarrier_h
#define G4He6GEMCoulombBarrier_h 1

#include "G4GEMCoulombBarrierHE.hh"
#include "globals.hh"

class G4He6GEMCoulombBarrier : public G4GEMCoulombBarrierHE
{
public:
  G4He6GEMCoulombBarrier() : G4GEMCoulombBarrierHE(6,2) {};
  ~G4He6GEMCoulombBarrier() {};

private:
  G4He6GEMCoulombBarrier(const G4He6GEMCoulombBarrier & right);

  const G4He6GEMCoulombBarrier & operator=(const G4He6GEMCoulombBarrier & right);
  G4bool operator==(const G4He6GEMCoulombBarrier & right) const;
  G4bool operator!=(const G4He6GEMCoulombBarrier & right) const;
  

};

#endif
