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
/// \file exoticphysics/phonon/include/XAluminumElectrodeSensitivity.hh
/// \brief Definition of the XAluminumElectrodeSensitivity class
//
//
// 20150818  Improve MT compatibility; hits collection should not be static

#ifndef XAluminumElectrodeSensitivity_h
#define XAluminumElectrodeSensitivity_h 1

#include "G4VSensitiveDetector.hh"
#include "XAluminumElectrodeHit.hh"
#include <iosfwd>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class XAluminumElectrodeSensitivity : public G4VSensitiveDetector {
public:
  XAluminumElectrodeSensitivity(const G4String&);
  virtual ~XAluminumElectrodeSensitivity();
  
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  
  XAluminumElectrodeHitsCollection* GetHitsCollection();

protected:
  void WriteHitInfo(const XAluminumElectrodeHit* aHit);

private:
  XAluminumElectrodeHitsCollection* fHitsCollection;

  static std::fstream* fWriter;         // For hit position output (temporary)
  static std::fstream* fWriter2;        // For hit timing/energy (temporary)
  
  G4int fHCID;                // Index of collection in event
};

#endif

