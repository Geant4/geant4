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
// $Id: PrimaryGeneratorAction.hh,v 1.1 2010-05-16 21:32:31 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
   ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

  public:
    
    G4ThreeVector GetNewUz()  {return newUz;};
    G4double GetAlphaMax()    {return alphaMax;};
        
    G4double RejectAccept();
    G4double InverseCumul();            
    
  private:
    G4ParticleGun*         particleGun;
 
    G4ThreeVector          newUz;   
    G4double               alphaMax;	//max angle
    
    G4int                  nPoints;	//tabulated function
    std::vector<G4double>  x;
    std::vector<G4double>  f;		//f(x)
    std::vector<G4double>  a;		//slopes
    std::vector<G4double>  Fc;		//cumulative of f
    G4double               fMax;	//max(f)

  private:
    void InitFunction();                        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
