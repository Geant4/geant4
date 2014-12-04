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
#ifndef Tst75HistoSet_h
#define Tst75HistoSet_h

#include "TstHistoSet.hh"

#include <vector>

#include "TProfile.h"

// class G4VParticleChange;

class Tst75HistoSet : public TstHistoSet
{

public:

   Tst75HistoSet( G4String, G4float );
   ~Tst75HistoSet();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:
   
   // data members
   //
   // plots for secondary proton
   //
   TH1F* fProton45;
   TH1F* fProton60;
   TH1F* fProton72;
   TH1F* fProton84;
   TH1F* fProton90;
   TH1F* fProton135;
   TH1F* fProKE;
   TH1F* fProCosTh;
   //
   // plots for secondary pi-
   //
   TH1F* fPim28deg;
   TH1F* fPim44deg;
   //
   // plots for secondary pi+
   //
   TH1F* fPip28deg;
   TH1F* fPip44deg;
   
   G4float fMaxKE;
   
   // G4VParticleChange* fInteraction;

};

#endif
