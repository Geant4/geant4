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
#ifndef Tst23ParticleChange_H
#define Tst23ParticleChange_H 1

#include "G4VParticleChange.hh"

class Tst23ParticleChange : public G4VParticleChange
{

   public:
   
      Tst23ParticleChange() : fIsFirstInter(false), fIncomingTrack(0) {}
      Tst23ParticleChange( bool isFirst ) : fIsFirstInter( isFirst ), fIncomingTrack(0) {}
      virtual ~Tst23ParticleChange() {}
      
      void           SetIncomingTrack( G4Track* trk ) { fIncomingTrack=trk; return; }
      bool           IsFisrtInteraction() const { return fIsFirstInter; }
      const G4Track* GetIncomingTrack()   const { return fIncomingTrack; }

   private:
   
      bool     fIsFirstInter;
      G4Track* fIncomingTrack;

};

#endif
