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
#ifndef Tst23ActionInit_h
#define Tst23ActionInit_h 1

#include "G4VUserActionInitialization.hh"

// forward declarations
//
class G4VUserPrimaryGeneratorAction;
class G4UserSteppingAction;

class Tst23ActionInit : public G4VUserActionInitialization
{

   public:
      
      // ctor & dtor
      Tst23ActionInit(); // : fPrimaryGen(0), fSteppingAct(0) {}
      virtual ~Tst23ActionInit() {}
      
      virtual void Build() const;
      virtual void BuildForMaster() const { return; }

      void SetAct( G4VUserPrimaryGeneratorAction* pg ){ fPrimaryGen=pg; return; }
      void SetAct( G4UserSteppingAction* sa )         { fSteppingAct=sa; return; }

   private:
            
      G4VUserPrimaryGeneratorAction* fPrimaryGen;
      G4UserSteppingAction*          fSteppingAct;
      
};

#endif
