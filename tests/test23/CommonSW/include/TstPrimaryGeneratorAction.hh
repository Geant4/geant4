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
#ifndef TstPrimaryGeneratorAction_h
#define TstPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4LorentzVector.hh"

class G4ParticleGun;
class G4Event;
class TstReader;
class Beam;


class TstPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction 
{

   public:

      // ctor & dtor
      TstPrimaryGeneratorAction() : fPartGun(0), fConfigPtr(0), fXSecOnTarget(0.), fIsInit(false) {} ;
      ~TstPrimaryGeneratorAction();  

      void  InitBeam( TstReader* );
      G4double               GetXSecOnTarget() const { return fXSecOnTarget; }
      const G4LorentzVector& GetLabV()         const { return fLabV; }
      const G4LorentzVector& GetLabP()         const { return fLabP; }
     
      void GeneratePrimaries( G4Event* );

   private:

      G4ParticleGun*  fPartGun;
      TstReader*      fConfigPtr;
      G4LorentzVector fLabV;
      G4LorentzVector fLabP;
      G4double        fXSecOnTarget;
      G4bool          fIsInit;

};

#endif
