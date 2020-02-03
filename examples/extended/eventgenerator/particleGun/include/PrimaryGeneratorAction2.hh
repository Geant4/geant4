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
/// \file eventgenerator/particleGun/include/PrimaryGeneratorAction2.hh
/// \brief Definition of the PrimaryGeneratorAction2 class
//
//
<<<<<<< HEAD
// $Id: PrimaryGeneratorAction2.hh 68024 2013-03-13 13:42:01Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef PrimaryGeneratorAction2_h
#define PrimaryGeneratorAction2_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction2
{
  public:
    PrimaryGeneratorAction2(G4ParticleGun*);    
   ~PrimaryGeneratorAction2();

  public:
   void GeneratePrimaries(G4Event*);

  public:        
    G4double RejectAccept();
    G4double InverseCumul();            
    
  private:    
    G4ParticleGun*         particleGun;
 
    G4int                  nPoints;     //tabulated function
    std::vector<G4double>  x;
    std::vector<G4double>  f;           //f(x)
    std::vector<G4double>  a;           //slopes
    std::vector<G4double>  Fc;          //cumulative of f
    G4double               fMax;        //max(f)

  private:
    void InitFunction();                        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
