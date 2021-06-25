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
/// \file eventgenerator/pythia/pythia8decayer/include/Py8DecayerPhysics
/// \brief Definition of the Py8DecayerPhysics class
///
/// \author J. Yarba; FNAL
///

#ifndef Py8DecayerPhysics_H
#define Py8DecayerPhysics_H

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Decay;

class Py8DecayerPhysics: public G4VPhysicsConstructor
{
  public:
    Py8DecayerPhysics(G4int verb=1);
    virtual ~Py8DecayerPhysics();

  protected:
    // methods
    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
