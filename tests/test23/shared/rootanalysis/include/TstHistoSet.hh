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
#ifndef TstHistoSet_h
#define TstHistoSet_h

#include "G4LorentzVector.hh"
#include "TstReader.hh"

// fwd declaration
class G4VParticleChange;
class G4Track;
class G4DecayProducts;

class TstHistoSet
{

   public:

//      enum NuERange        { kNone=-1, kR_0_2=0, kR_2_5=1, kR_5_10=2, kR_10_20=3, kR_20_50=4 };
//      enum InteractionType { kOther=-1, kFromP=0, kFromPim=1, kFromPip=2 };
   
      TstHistoSet(G4String) : fDoResDecay(false) {}
      virtual ~TstHistoSet() {}
      
      virtual void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& ) = 0;
      virtual void Write( G4int, G4double ) = 0;
      
      void SetDoResDecay( G4bool rdec ) { fDoResDecay=rdec; return; }
      
   protected:
   
      void AccountForResDecay( G4VParticleChange* );
//      void AccountForPionDecay( const G4Track* );

      G4bool           fDoResDecay;
//      NuERange         fNuERange;
//      NuERange         fNubarERange;
//      G4DecayProducts* fPionDecay;

};

#endif
