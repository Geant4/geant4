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
// GEANT4 tag $Name:  $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4He3Builder
//
// Author: 2013 P. Arce
//
// Modified
//// 12.04.2017 A.Dotti move to new design with base class
//----------------------------------------------------------------------------
//
#ifndef G4He3Builder_h
#define G4He3Builder_h 1

#include "G4PhysicsBuilderInterface.hh"
#include "globals.hh"

#include "G4He3InelasticProcess.hh"
#include "G4VHe3Builder.hh"
#include <vector>

class G4He3Builder : public G4PhysicsBuilderInterface
{
  public: 
    G4He3Builder();
    virtual ~G4He3Builder() {}

    virtual void Build() final override;
    virtual void RegisterMe(G4PhysicsBuilderInterface * aB) final override;

  private:
    G4He3InelasticProcess * theHe3Inelastic;
    
    std::vector<G4VHe3Builder *> theModelCollections;

    G4bool wasActivated;
};

#endif

