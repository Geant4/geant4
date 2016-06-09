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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraPrimaryGeneratorAction.hh
//   **********************************************
//  
//    Class used in the definition of the optical photons source
//    A plane, circular source is used. Depending on the source position, optical
//    photons may reach the UVscope directly or after reflection. By default direct
//    incidence is used. The source parameters can be set directly in this class
//    or through the GeneralParticleSource  messenger class.
//
#ifndef UltraPrimaryGeneratorAction_H
#define UltraPrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

class G4GeneralParticleSource;
class G4Event;

class UltraPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    UltraPrimaryGeneratorAction();
    ~UltraPrimaryGeneratorAction();

  public:

    void GeneratePrimaries(G4Event* anEvent);

    G4GeneralParticleSource* GetParticleGun(){return particleGun;}    

  private:

    G4GeneralParticleSource* particleGun;

};

#endif 


