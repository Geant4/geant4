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
/// \file electromagnetic/TestEm1/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
<<<<<<< HEAD
// $Id: DetectorConstruction.hh 84815 2014-10-21 12:19:02Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
<<<<<<< HEAD
   ~DetectorConstruction();

  public:
=======
    virtual ~DetectorConstruction();
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
     
    void SetSize     (G4double);              
    void SetMaterial (const G4String&);
  
    inline const G4VPhysicalVolume* GetWorld() const {return fPBox;};
    inline G4double GetSize() const                  {return fBoxSize;};      
    inline const G4Material* GetMaterial() const     {return fMaterial;};
     
<<<<<<< HEAD
     void               PrintParameters();
                       
  private:
  
     G4VPhysicalVolume*    fPBox;
     G4LogicalVolume*      fLBox;
=======
    void   PrintParameters();
    void   DefineMaterials();
                       
  private:
  
    G4VPhysicalVolume*    fPBox;
    G4LogicalVolume*      fLBox;
    G4Box*                fBox;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
     
    G4double              fBoxSize;
    G4Material*           fMaterial;

<<<<<<< HEAD
  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
=======
    DetectorMessenger* fDetectorMessenger;
    G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

