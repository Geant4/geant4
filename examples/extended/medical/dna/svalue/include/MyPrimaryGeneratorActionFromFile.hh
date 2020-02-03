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
// Derived from 
//  https://twiki.cern.ch/twiki/bin/view/Geant4/QuickMigrationGuideForGeant4V10
// Courtesy of A. Dotti
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/include/MyPrimaryGeneratorActionFromFile.hh
/// \brief Declaration of the MyPrimaryGeneratorActionFromFile class

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4VStateDependent.hh"

class MyFileReader;
class DetectorConstruction;

class MyPrimaryGeneratorActionFromFile
 : public G4VUserPrimaryGeneratorAction,
   public G4VStateDependent
{
 public:
   MyPrimaryGeneratorActionFromFile();
   virtual ~MyPrimaryGeneratorActionFromFile();
   
  virtual G4bool Notify(G4ApplicationState requestedState);
  
  virtual void GeneratePrimaries(G4Event* anEvent);
 
  G4ParticleGun* GetParticleGun() const
  {
    return fParticleGun;
  }

 private:
   static MyFileReader* fileReader;
   G4ParticleGun* fParticleGun;
   const DetectorConstruction* fDetector;
};
