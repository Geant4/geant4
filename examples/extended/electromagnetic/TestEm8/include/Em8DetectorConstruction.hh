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
// $Id: Em8DetectorConstruction.hh,v 1.7 2004/05/27 08:39:05 grichine Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

#ifndef Em8DetectorConstruction_h
#define Em8DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class G4UniformMagField;
class Em8DetectorMessenger;
class Em8CalorimeterSD;



class Em8DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em8DetectorConstruction();
   ~Em8DetectorConstruction();

  public:
     
  void SetAbsorberMaterial (G4String);     
  void SetAbsorberThickness(G4double);     
  void SetAbsorberRadius(G4double);          
      
  void SetAbsorberZpos(G4double);

  void SetWorldMaterial(G4String);
  void SetWorldSizeZ(G4double);
  void SetWorldSizeR(G4double);

  void SetGammaCut(G4double    cut){fGammaCut    = cut;};
  void SetElectronCut(G4double cut){fElectronCut = cut;};
  void SetPositronCut(G4double cut){fPositronCut = cut;};


  //  void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
  void PrintCalorParameters(); 
                    
  G4Material* GetWorldMaterial()    {return fWorldMaterial;};
  G4double GetWorldSizeZ()          {return fWorldSizeZ;}; 
  G4double GetWorldSizeR()          {return fWorldSizeR;};
     
  G4double GetAbsorberZ()        {return fAbsorberZ;}; 
  G4double GetStartAbsZ()           {return fStartAbsZ;};
  G4double GetEndAbsZ()             {return fEndAbsZ;};

  G4Material* GetAbsorberMaterial()  {return fAbsorberMaterial;};
  G4double    GetAbsorberThickness() {return fAbsorberThickness;};      
  G4double GetAbsorberRadius()       {return fAbsorberRadius;};
     
  const G4VPhysicalVolume* GetphysiWorld() {return fPhysicsWorld;};           
  const G4VPhysicalVolume* GetAbsorber()   {return fPhysicsAbsorber;};
  G4LogicalVolume* GetLogicalAbsorber()    {return fLogicAbsorber;};
                 
  private:
     
  static const G4double fDelta;

  G4bool             fWorldChanged;
  G4double           fAbsorberThickness;
  G4double           fAbsorberRadius;

  G4Material*        fWindowMat ;
  G4double           fWindowThick ;


 
  G4double           fAbsorberZ ;
  G4double           fStartAbsZ , fEndAbsZ ;
     
  G4Material*        fWorldMaterial;
  G4double           fWorldSizeR;
  G4double           fWorldSizeZ;
            
  G4Tubs*            fSolidWorld;     //pointer to the solid World 
  G4LogicalVolume*   fLogicWorld;     //pointer to the logical World
  G4VPhysicalVolume* fPhysicsWorld;   //pointer to the physical World


  G4Material*        fAbsorberMaterial;
  G4Tubs*            fSolidAbsorber;    //pointer to the solid Absorber
  G4LogicalVolume*   fLogicAbsorber;   //pointer to the logical Absorber
  G4VPhysicalVolume* fPhysicsAbsorber; //pointer to the physical Absorber

  G4double fElectronCut, fGammaCut, fPositronCut;
     
  Em8DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
  Em8CalorimeterSD*     fCalorimeterSD;      //pointer to the sensitive detector
  G4Region*             fRegGasDet;
      
  private:
    
  void DefineMaterials();
  void ComputeCalorParameters();
  G4VPhysicalVolume* ConstructCalorimeter();
     
};

//////////////////////////////////////////////////////////////////////

inline void Em8DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter

  if( !fWorldChanged )
  {
       //  WorldSizeR=2.*AbsorberRadius ;
       //  WorldSizeZ=2.*AbsorberThickness ;
  } 
  fWorldSizeZ = fAbsorberThickness + 2*fWindowThick + 2*fDelta;
  fWorldSizeR = fAbsorberRadius + fDelta;

  fStartAbsZ = fAbsorberZ - 0.5*fAbsorberThickness; 
  fEndAbsZ   = fAbsorberZ + 0.5*fAbsorberThickness; 

}

#endif

