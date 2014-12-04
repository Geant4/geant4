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
#ifndef ExecProcessLevel_H 
#define ExecProcessLevel_H 1

#include "ExecBase.hh"
#include "Beam.hh"
#include "TstTarget.hh"

// fwd declarations
// class TstTarget;
class G4Region;
class G4HadronicProcess;
class G4ProcessManager;
class G4Track;
class G4Step;
class G4StepPoint;

class ExecProcessLevel : public ExecBase
{

   public:
         
      ExecProcessLevel( const TstReader* pset ) : ExecBase(pset), fXSecOnTarget(0.) {}
      virtual ~ExecProcessLevel();
      
      virtual G4VParticleChange* DoEvent();
      
      const TstTarget* GetTarget() const { return fTarget; } 
      const Beam*      GetBeam()   const { return fBeam; }
      G4double         GetXSecOnTarget() const { return fXSecOnTarget; }
      
   protected:
   
      ExecProcessLevel();
      virtual void InitSetup(   const TstReader* );
      virtual void InitBeam(    const TstReader* );
      virtual void InitProcess( const TstReader* ) = 0;

      // protected data member
      G4HadronicProcess* fProcWrapper;
      G4double           fXSecOnTarget;
   
   private:
         
      // data members
      TstTarget*         fTarget;
      G4Region*          fRegion;
      G4ProcessManager*  fProcManager;
      Beam*              fBeam;
      G4Track*           fTrack;
      G4Step*            fStep;
      G4StepPoint*       fStepPoint;
      G4VParticleChange* fPartChange;

};

#endif
