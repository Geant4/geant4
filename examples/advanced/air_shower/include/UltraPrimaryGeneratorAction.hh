//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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


