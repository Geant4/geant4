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
/// \file InterPulseAction.hh
/// \brief Definition of the InterPulseAction class

// author: Le Tuan Anh

#ifndef InterPulseAction_h
#define InterPulseAction_h 1

#include "PulseAction.hh"
#include "globals.hh"

#include <vector>
class G4Track;

class InterPulseAction: public PulseAction
{
public:
  explicit InterPulseAction(const G4String& pulse, G4bool useHisto=false,
                  G4double pulsePeriod = 0, G4int npulses = 0);
  ~InterPulseAction() override = default;
  void PreUserTrackingAction(const G4Track *) override;
  void SetPulsePeriod(G4double tp) {fPulsePeriod = tp;}
  G4double GetPulsePeriod() const{return fPulsePeriod;}
  G4int GetNumberOfPulse() const{return fNumberOfPulse;}
private:
    G4int WhichPulse() const;
    G4int fNumberOfPulse = 1;//this must be always = 1 if no pulse
    G4double fPulsePeriod = 0;
};

#endif

