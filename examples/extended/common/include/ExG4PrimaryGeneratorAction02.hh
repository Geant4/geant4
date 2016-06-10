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
// $Id$
//
/// \file ExG4PrimaryGeneratorAction02.hh
/// \brief Definition of the ExG4PrimaryGeneratorAction02 class
//

#ifndef ExG4PrimaryGeneratorAction02_h
#define ExG4PrimaryGeneratorAction02_h 1

/// \file ExG4PrimaryGeneratorAction02.hh
/// \brief Definition of the ExG4PrimaryGeneratorAction02 class 

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

/// \ingroup primary_generator
/// \brief The primary generator class with general particle source
///
/// \author I. Hrivnacova; IPN Orsay

class ExG4PrimaryGeneratorAction02 : public G4VUserPrimaryGeneratorAction
{
  public:
    ExG4PrimaryGeneratorAction02();    
    ~ExG4PrimaryGeneratorAction02();

    // methods
    virtual void GeneratePrimaries(G4Event*);

  private:
    // static data members
    static const G4String fgkDefaultParticleName;
    static const G4double fgkDefaultEnergy;

    // data members
    G4GeneralParticleSource*  fGeneralParticleSource; //pointer a to G4 service class
};

#endif


