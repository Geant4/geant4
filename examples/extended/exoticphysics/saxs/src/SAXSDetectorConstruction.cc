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
/// \file SAXSDetectorConstruction.cc
/// \brief Implementation of the SAXSDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SAXSDetectorConstruction.hh"
#include "SAXSSensitiveDetector.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4IntersectionSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UniformMagField.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDParticleFilter.hh"

#include <cmath>
#include <sstream>
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4ExtendedMaterial.hh"
#include "G4MIData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSDetectorConstruction::SAXSDetectorConstruction():
  G4VUserDetectorConstruction(), fWorldLogic(0)
{          
  G4cout << "### DetectorConstruction Instantiated ###" << G4endl;

  //instantiate the messenger (set methods will be called after construct)
  fMessenger = new SAXSDetectorConstructionMessenger(this);
        
  //set geometrical variables
  SetGeometricalVariables();        
    
  //Initialization           
  fPhantomMaterialIndex = 1;   
    
  fComp0 = 0.0;    //components of the "Medical Material (MedMat)"
  fComp1 = 1.0;    //Through macro I can set one Medical Material only
  fComp2 = 0.0;
  fComp3 = 0.0;
    
  fCustomMatDensity = 1.00; //g/mol
  fCustomMatHmassfract = 0.1119;
  fCustomMatCmassfract = 0.;
  fCustomMatNmassfract = 0.;
  fCustomMatOmassfract = 0.8881;
  fCustomMatNamassfract = 0.;
  fCustomMatPmassfract = 0.;
  fCustomMatSmassfract = 0.;
  fCustomMatClmassfract = 0.;
  fCustomMatKmassfract = 0.;
  fCustomMatCamassfract = 0.;
    
  fCustomMatFF = "";  //MIFF filename for a custom (extended) material
    
  fSensitiveVolume = 0;
    
  fIWantSlits = false;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

SAXSDetectorConstruction::~SAXSDetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSDetectorConstruction::DefineMaterials()
{
  //Define the NIST manager
  G4NistManager* NistMan = G4NistManager::Instance();
    
  //Define the required elements for compounds
  G4Element* elH  = NistMan->FindOrBuildElement("H");
  G4Element* elC  = NistMan->FindOrBuildElement("C");
  G4Element* elN  = NistMan->FindOrBuildElement("N");
  G4Element* elO  = NistMan->FindOrBuildElement("O");
  G4Element* elNa = NistMan->FindOrBuildElement("Na"); 
  G4Element* elP  = NistMan->FindOrBuildElement("P");  
  G4Element* elS  = NistMan->FindOrBuildElement("S");  
  G4Element* elCl = NistMan->FindOrBuildElement("Cl");  
  G4Element* elK = NistMan->FindOrBuildElement("K");  
  G4Element* elCa = NistMan->FindOrBuildElement("Ca");      
        
  //variable definition
  G4double d;                          //density
  G4int nel;                           //number of elements
  G4String matname;
    
  //Air
  d = 1.29*mg/cm3;
  nel = 2;
  G4double tAir = 293.15 * CLHEP::kelvin;  //20° Celsius
  G4double pAir = 1.*CLHEP::atmosphere;          //1 atm
  fAir = new G4Material("Air", d, nel, kStateGas, tAir, pAir);
  fAir->AddElement(elN, 0.7);
  fAir->AddElement(elO, 0.3);
                  
  //Fat (Tartari2002) (FF from Tartari2002)
  G4double d_Fat = 0.923*g/cm3;
  nel = 3;
  matname = "Fat_MI";
  fFat = new G4Material(matname, d_Fat, nel);
  fFat->AddElement(elH, 0.119);
  fFat->AddElement(elC, 0.772);
  fFat->AddElement(elO, 0.109);
    
  //Water (FF from Tartari2002)
  G4double d_Water = 1.*g/cm3;
  nel = 2;
  matname = "Water_MI";
  fWater = new G4Material(matname, d_Water, nel);
  fWater->AddElement(elH, 2);
  fWater->AddElement(elO, 1);
    
  //BoneMatrix (Collagen) (FF from Tartari2002)
  G4double d_BoneMatrix = 1.263*g/cm3; 
  nel = 4;
  matname = "BoneMatrix_MI";
  fBoneMatrix = new G4Material(matname, d_BoneMatrix, nel);
  fBoneMatrix->AddElement(elH, 0.0344);
  fBoneMatrix->AddElement(elC, 0.7140);
  fBoneMatrix->AddElement(elN, 0.1827);
  fBoneMatrix->AddElement(elO, 0.0689);
    
  //Mineral (Hydroxyapatite) (Tartari2002) (FF from Tartari2002)
  G4double d_Mineral = 2.74*g/cm3;
  nel = 4;
  matname = "Mineral_MI";
  fMineral = new G4Material(matname, d_Mineral, nel);
  fMineral->AddElement(elH, 0.002);
  fMineral->AddElement(elO, 0.414);
  fMineral->AddElement(elP, 0.185);
  fMineral->AddElement(elCa, 0.399);
    
  //Medical Material (compostion of Water, Fat, BoneMatrix and Mineral)
  G4double comp[] = {fComp0, fComp1, fComp2, fComp3};
  G4double d_MedMat =
    1/(comp[0]/d_Fat+comp[1]/d_Water+comp[2]/d_BoneMatrix+comp[3]/d_Mineral);
  G4int n_MedMat = 0;
  for (size_t i=0; i<4; i++) {
    if (comp[i]>0) n_MedMat++;
    if (comp[i]<0 || comp[i]>1) {
      G4String excep =
        "Error in Medical Material composition: comp[i]<0 or comp[i]>1";
      G4Exception("DetectorConstuction::DefineMaterials()","dc0001",
                  FatalException,excep);
      return;
    }
  }        
  std::stringstream ss0,ss1,ss2,ss3;
  ss0 << comp[0];        
  ss1 << comp[1];        
  ss2 << comp[2];        
  ss3 << comp[3];        
  if (comp[0]==0 || comp[0]==1) ss0 << ".00";
  if (comp[1]==0 || comp[1]==1) ss1 << ".00";
  if (comp[2]==0 || comp[2]==1) ss2 << ".00";
  if (comp[3]==0 || comp[3]==1) ss3 << ".00";
  if (ss0.str().size()<4) ss0 << "0";
  if (ss1.str().size()<4) ss1 << "0";
  if (ss2.str().size()<4) ss2 << "0";
  if (ss3.str().size()<4) ss3 << "0";
  if (ss0.str().size()!=4 || ss1.str().size()!=4 || ss2.str().size()!=4 
      || ss3.str().size()!=4) {
    G4String excep = 
      "Error in MedMaterial composition: check the digits of the elements of comp";
    G4Exception("DetectorConstuction::DefineMaterials()","dc0002",
                FatalException,excep);
    return;
  }
  matname = "MedMat_"+ss0.str()+"_"+ss1.str()+"_"+ss2.str()+"_"+ss3.str();
  fMedMat = new G4Material(matname, d_MedMat, n_MedMat);
  if (comp[0]) fMedMat->AddMaterial(fFat, comp[0]);
  if (comp[1]) fMedMat->AddMaterial(fWater, comp[1]);
  if (comp[2]) fMedMat->AddMaterial(fBoneMatrix, comp[2]);
  if (comp[3]) fMedMat->AddMaterial(fMineral, comp[3]);        
  if (comp[0]+comp[1]+comp[2]+comp[3] != 1) {
    G4String excep = "Error in Medical Material composition: sum(comp) != 1";
    G4Exception("DetectorConstuction::DefineMaterials()",
                "dc0003",FatalException,excep);
    return;
  }
  //If the user wants to use more than one MedMat, he has to create the mix 
  //and label the material properly, such as "MedMat_0.55_0.25_0.05_0.15". 
  //Such a name enables the automatic form factor calculation.
    
  //PMMA (FF from Tartari2002)
  d = 1.18*g/cm3;
  nel = 3;
  matname = "PMMA_MI";
  fPMMA = new G4Material(matname, d, nel);
  fPMMA->AddElement(elH, 8);
  fPMMA->AddElement(elC, 5);
  fPMMA->AddElement(elO, 2);
    
  //Adipose (Poletti2002) (FF from Poletti2002)
  d = 0.92*g/cm3;
  nel = 4;
  matname = "adipose_MI";
  fAdipose = new G4Material(matname, d, nel);
  fAdipose->AddElement(elH, 0.124);
  fAdipose->AddElement(elC, 0.765);
  fAdipose->AddElement(elN, 0.004);
  fAdipose->AddElement(elO, 0.107);
    
  //Glandular (Poletti2002) (FF from Poletti2002)
  d = 1.04*g/cm3;
  nel = 4;
  matname = "glandular_MI";
  fGlandular = new G4Material(matname, d, nel);
  fGlandular->AddElement(elH, 0.093);
  fGlandular->AddElement(elC, 0.184);
  fGlandular->AddElement(elN, 0.044);
  fGlandular->AddElement(elO, 0.679);
    
  //human breast 50/50 (ICRU44) (FF from Peplow1998)
  d = 0.96*g/cm3;
  nel = 3;
  matname = "breast5050_MI";
  fBreast5050 = new G4Material(matname, d, nel);
  fBreast5050->AddElement(elH, 0.115);
  fBreast5050->AddElement(elC, 0.387);
  fBreast5050->AddElement(elO, 0.498);
    
  //Liver (ICRU46) (FF_pork_liver_Peplow1998)
  d = 1.06*g/cm3; 
  nel = 9;
  matname = "liver_MI";
  fliver = new G4Material(matname, d, nel);
  fliver->AddElement(elH, 0.102);
  fliver->AddElement(elC, 0.139);
  fliver->AddElement(elN, 0.030);
  fliver->AddElement(elO, 0.716);
  fliver->AddElement(elNa, 0.002);
  fliver->AddElement(elP, 0.003);
  fliver->AddElement(elS, 0.003);
  fliver->AddElement(elCl, 0.002);
  fliver->AddElement(elK, 0.003);
    
  //Kidney (ICRU46) (FF_pork_kidney_Peplow1998)
  d = 1.05*g/cm3; 
  nel = 10;
  matname = "kidney_MI";
  fkidney = new G4Material(matname, d, nel);
  fkidney->AddElement(elH, 0.103);
  fkidney->AddElement(elC, 0.132);
  fkidney->AddElement(elN, 0.030);
  fkidney->AddElement(elO, 0.724);
  fkidney->AddElement(elNa, 0.002);
  fkidney->AddElement(elP, 0.002);
  fkidney->AddElement(elS, 0.002);
  fkidney->AddElement(elCl, 0.002);
  fkidney->AddElement(elK, 0.002);
  fkidney->AddElement(elCa, 0.001);  
    
  //Lexan (Polycarbonate) (FF from Peplow1998)
  d = 1.221*g/cm3;
  nel = 3;
  matname = "Lexan_MI";
  fLexan = new G4Material(matname, d, nel);
  fLexan->AddElement(elH, 14);
  fLexan->AddElement(elC, 16);
  fLexan->AddElement(elO, 3);
    
  //Kapton (FF from Peplow1998)
  d = 1.42*g/cm3;
  nel = 4;
  matname = "Kapton_MI";
  fKapton = new G4Material(matname, d, nel);
  fKapton->AddElement(elH, 28);
  fKapton->AddElement(elC, 35);
  fKapton->AddElement(elN, 2);
  fKapton->AddElement(elO, 7);
    
  //Carcinoma (muscle ICRU44) (FF from Kidane1999)
  d = 1.05*g/cm3; //check the density
  nel = 9;
  matname = "carcinoma_MI";
  fcarcinoma = new G4Material(matname, d, nel);
  fcarcinoma->AddElement(elH, 0.102);
  fcarcinoma->AddElement(elC, 0.143);
  fcarcinoma->AddElement(elN, 0.034);
  fcarcinoma->AddElement(elO, 0.710);
  fcarcinoma->AddElement(elNa, 0.001);
  fcarcinoma->AddElement(elP, 0.002);
  fcarcinoma->AddElement(elS, 0.003);
  fcarcinoma->AddElement(elCl, 0.001);
  fcarcinoma->AddElement(elK, 0.004);
    
  //Nylon (FF from Kosanetzky1987)
  d = 1.15*g/cm3;
  nel = 4;
  matname = "Nylon_MI";
  fNylon = new G4Material(matname, d, nel);          
  fNylon->AddElement(elH, 11);
  fNylon->AddElement(elC, 6);
  fNylon->AddElement(elN, 1);
  fNylon->AddElement(elO, 1);
    
  //Polyethylene (FF from Kosanetzky1987) 
  d = 0.94*g/cm3;        //MDPE => 0.92*g/cm3, HDPE => 0.94*g/cm3
  nel = 2;
  matname = "Polyethylene_MI";
  fPolyethylene = new G4Material(matname, d, nel);
  fPolyethylene->AddElement(elH, 4);
  fPolyethylene->AddElement(elC, 2);
    
  //Polystyrene (FF from Kosanetzky1987)
  d = 1.05*g/cm3;
  nel = 2;
  matname = "Polystyrene_MI";
  fPolystyrene = new G4Material(matname, d, nel);
  fPolystyrene->AddElement(elH, 8);
  fPolystyrene->AddElement(elC, 8);          
    
  //GrayMatter (DeFelici2008) (FF from DeFelici2008)
  d = 0.991*g/cm3; 
  nel = 3;
  matname = "grayMatter_MI";
  fGrayMatter = new G4Material(matname, d, nel);
  fGrayMatter->AddElement(elH, 0.1127);
  fGrayMatter->AddElement(elC, 0.0849);
  fGrayMatter->AddElement(elO, 0.8024);
    
  //WhiteMatter (DeFelici2008) (FF from DeFelici2008)
  d = 0.983*g/cm3; 
  nel = 3;
  matname = "whiteMatter_MI";
  fWhiteMatter = new G4Material(matname, d, nel);
  fWhiteMatter->AddElement(elH, 0.1134);
  fWhiteMatter->AddElement(elC, 0.1621);
  fWhiteMatter->AddElement(elO, 0.7245);
    
  //Blood (beef) (ICRU46) (FF from Peplow1998)
  d = 1.06*g/cm3; 
  nel = 9;
  matname = "blood_MI";
  fbeefBlood = new G4Material(matname, d, nel);
  fbeefBlood->AddElement(elH, 0.102);
  fbeefBlood->AddElement(elC, 0.11);
  fbeefBlood->AddElement(elN, 0.033);
  fbeefBlood->AddElement(elO, 0.746);
  fbeefBlood->AddElement(elNa, 0.001);
  fbeefBlood->AddElement(elP, 0.001);
  fbeefBlood->AddElement(elS, 0.002);
  fbeefBlood->AddElement(elCl, 0.003);
  fbeefBlood->AddElement(elK, 0.002);
    
  //Formaline (FF from Peplow1998)
  d = 1.083*g/cm3; 
  nel = 3;
  matname = "Formaline_MI";
  fFormaline = new G4Material(matname, d, nel);
  fFormaline->AddElement(elH, 2);
  fFormaline->AddElement(elC, 1);
  fFormaline->AddElement(elO, 1);          
    
  //Acetone (FF from Cozzini2010)
  d = 0.7845*g/cm3; 
  nel = 3;
  matname = "Acetone_MI";
  fAcetone = new G4Material(matname, d, nel);
  fAcetone->AddElement(elH, 6);
  fAcetone->AddElement(elC, 3);
  fAcetone->AddElement(elO, 1);
    
  //Hperoxide (FF from Cozzini2010)
  d = 1.11*g/cm3; 
  nel = 2;
  matname = "Hperoxide_MI";
  fHperoxide = new G4Material(matname, d, nel);
  fHperoxide->AddElement(elH, 2);
  fHperoxide->AddElement(elO, 2);
    
  //CIRS30-70 (Poletti2002) (FF from Poletti2002)
  d = 0.97*g/cm3; 
  nel = 5;
  matname = "CIRS30-70_MI";
  fCIRS3070 = new G4Material(matname, d, nel);
  fCIRS3070->AddElement(elH, 0.1178);
  fCIRS3070->AddElement(elC, 0.7512);
  fCIRS3070->AddElement(elN, 0.0066);
  fCIRS3070->AddElement(elO, 0.1214);
  fCIRS3070->AddElement(elCa, 0.0030);
    
  //CIRS50-50 (Poletti2002) (FF from Poletti2002)
  d = 0.98*g/cm3; 
  nel = 5;
  matname = "CIRS50-50_MI";
  fCIRS5050 = new G4Material(matname, d, nel);
  fCIRS5050->AddElement(elH, 0.1110);
  fCIRS5050->AddElement(elC, 0.7274);
  fCIRS5050->AddElement(elN, 0.0104);
  fCIRS5050->AddElement(elO, 0.1482);
  fCIRS5050->AddElement(elCa, 0.0030);
    
  //CIRS70-30 (Poletti2002) (FF from Poletti2002)
  d = 1.01*g/cm3; 
  nel = 5;
  matname = "CIRS70-30_MI";
  fCIRS7030 = new G4Material(matname, d, nel);
  fCIRS7030->AddElement(elH, 0.1172);
  fCIRS7030->AddElement(elC, 0.7378);
  fCIRS7030->AddElement(elN, 0.0130);
  fCIRS7030->AddElement(elO, 0.1244);
  fCIRS7030->AddElement(elCa, 0.0076);
          
  //RMI454 (Poletti2002) (FF from Poletti2002)
  d = 0.98*g/cm3; 
  nel = 4;
  matname = "RMI454_MI";
  fRMI454 = new G4Material(matname, d, nel);
  fRMI454->AddElement(elH, 0.0924);
  fRMI454->AddElement(elC, 0.6935);
  fRMI454->AddElement(elN, 0.0198);
  fRMI454->AddElement(elO, 0.1943);
        
  //Bone (King2011 decomposition) (FF from King2011) 
  d = 1.344*g/cm3; 
  nel = 6;
  matname = "bone_MI";
  fBone = new G4Material(matname, d, nel);
  fBone->AddElement(elH, 0.0582);
  fBone->AddElement(elC, 0.3055);
  fBone->AddElement(elN, 0.0347);
  fBone->AddElement(elO, 0.3856);
  fBone->AddElement(elP, 0.0684);
  fBone->AddElement(elCa, 0.1476);

  //FatLowX (Tartari2002) (FF_fat_Tartari2002_joint_lowXdata_ESRF2003)
  nel = 3;
  matname = "FatLowX_MI";
  ffatLowX = new G4Material(matname, d_Fat, nel);
  ffatLowX->AddElement(elH, 0.119);
  ffatLowX->AddElement(elC, 0.772);
  ffatLowX->AddElement(elO, 0.109);
          
  //BonematrixLowX (Collagen)
  //(Tartari2002) (FF_bonematrix_Tartari2002_joint_lowXdata)
  nel = 4;
  matname = "BoneMatrixLowX_MI";
  fbonematrixLowX = new G4Material(matname, d_BoneMatrix, nel);
  fbonematrixLowX->AddElement(elH, 0.0344);
  fbonematrixLowX->AddElement(elC, 0.7140);
  fbonematrixLowX->AddElement(elN, 0.1827);
  fbonematrixLowX->AddElement(elO, 0.0689);
          
  //dryBoneLowX
  //(Tartari2002) (FF_dryBone_Tartari2002_joint_lowXdata_ESRF2003)
  d = 2.06*g/cm3; 
  nel = 6;
  matname = "dryBoneLowX_MI";
  fdryBoneLowX = new G4Material(matname, d, nel);
  fdryBoneLowX->AddElement(elH, 0.0112);
  fdryBoneLowX->AddElement(elC, 0.2013);
  fdryBoneLowX->AddElement(elN, 0.0515);
  fdryBoneLowX->AddElement(elO, 0.3148);
  fdryBoneLowX->AddElement(elP, 0.1327);
  fdryBoneLowX->AddElement(elCa, 0.2885);

  //CustomMat (FF read from file)
  nel = 0;
  G4double sumMF = 0.;
  G4cout << "CustomMat composition: " << G4endl;
  if (fCustomMatHmassfract) {
    G4cout << "CustomMatHmassfract: " << fCustomMatHmassfract << G4endl;
    nel++;
    sumMF += fCustomMatHmassfract;
  }
  if (fCustomMatCmassfract) { 
    G4cout << "CustomMatCmassfract: " << fCustomMatCmassfract << G4endl;
    nel++;
    sumMF += fCustomMatCmassfract;
  }
  if (fCustomMatNmassfract) {
    G4cout << "CustomMatNmassfract: " << fCustomMatNmassfract << G4endl;
    nel++;
    sumMF += fCustomMatNmassfract;
  }
  if (fCustomMatOmassfract) {
    G4cout << "CustomMatOmassfract: " << fCustomMatOmassfract << G4endl;
    nel++;
    sumMF += fCustomMatOmassfract;
  }
  if (fCustomMatNamassfract) {
    G4cout << "CustomMatNamassfract: " << fCustomMatNamassfract << G4endl;
    nel++;
    sumMF += fCustomMatNamassfract;
  }
  if (fCustomMatPmassfract) {
    G4cout << "CustomMatPmassfract: " << fCustomMatPmassfract << G4endl;
    nel++;
    sumMF += fCustomMatPmassfract;
  }
  if (fCustomMatSmassfract) {
    G4cout << "CustomMatSmassfract: " << fCustomMatSmassfract << G4endl;
    nel++;
    sumMF += fCustomMatSmassfract;
  }
  if (fCustomMatClmassfract) {
    G4cout << "CustomMatClmassfract: " << fCustomMatClmassfract << G4endl;
    nel++;
    sumMF += fCustomMatClmassfract;
  } 
  if (fCustomMatKmassfract) {
    G4cout << "CustomMatKmassfract: " << fCustomMatKmassfract << G4endl;
    nel++;
    sumMF += fCustomMatKmassfract;
  } 
  if (fCustomMatCamassfract) {
    G4cout << "CustomMatCamassfract: " << fCustomMatCamassfract << G4endl;
    nel++;
    sumMF += fCustomMatCamassfract;
  }
  if (sumMF == 0.) {
    //set a default material (water), 
    //otherwiswe an error appears in the interactive mode
    fCustomMatDensity = 1.00; //g/cm3
    fCustomMatHmassfract = 0.1119;
    G4cout << "CustomMat set, but not used!" << G4endl;
    G4cout << "CustomMatHmassfract: " << fCustomMatHmassfract << G4endl;
    fCustomMatOmassfract = 0.8881;
    G4cout << "CustomMatOmassfract: " << fCustomMatOmassfract << G4endl;
    nel = 2;
    sumMF = 1;                
  }        
  if (sumMF != 1.) {
    G4String excep =
      "Error in Custom Material composition: check elemental mass fractions";
    G4Exception("DetectorConstuction::DefineMaterials()",
                "dc0004",FatalException,excep);
    return;
  }   
        
  d = fCustomMatDensity*g/cm3; 
  matname = "CustomMat";
  //Notice: this is an extended material
  fCustomMat = new G4ExtendedMaterial(matname, d, nel);
  if (fCustomMatHmassfract)
    fCustomMat->AddElement(elH, fCustomMatHmassfract);        
  if (fCustomMatCmassfract)
    fCustomMat->AddElement(elC, fCustomMatCmassfract);
  if (fCustomMatNmassfract)
    fCustomMat->AddElement(elN, fCustomMatNmassfract);
  if (fCustomMatOmassfract)
    fCustomMat->AddElement(elO, fCustomMatOmassfract);
  if (fCustomMatNamassfract)
    fCustomMat->AddElement(elNa, fCustomMatNamassfract);
  if (fCustomMatPmassfract)
    fCustomMat->AddElement(elP, fCustomMatPmassfract);
  if (fCustomMatSmassfract)
    fCustomMat->AddElement(elS, fCustomMatSmassfract);
  if (fCustomMatClmassfract)
    fCustomMat->AddElement(elCl, fCustomMatClmassfract);
  if (fCustomMatKmassfract)
    fCustomMat->AddElement(elK, fCustomMatKmassfract);
  if (fCustomMatCamassfract)
    fCustomMat->AddElement(elCa, fCustomMatCamassfract);
  //Register MI extension        
  fCustomMat->RegisterExtension
    (std::unique_ptr<G4MIData>(new G4MIData("MI")));
  G4MIData* dataMICustomMat = (G4MIData*)fCustomMat->RetrieveExtension("MI");
  dataMICustomMat->SetFilenameFF(fCustomMatFF);

  //Nist Materials
  fLead = NistMan->FindOrBuildMaterial("G4_Pb"); 
  fTungsten = NistMan->FindOrBuildMaterial("G4_W");
  fGe = NistMan->FindOrBuildMaterial("G4_Ge");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSDetectorConstruction::SetGeometricalVariables()
{
  //World
  fWorldSize = 10000.*mm; 
    
  //Phantom
  fPhantomDiameter = 50.*mm;
  fPhantomHeight = 15.*mm;
  fPhantomZ = 0.*mm;
            
  //setup angle (rad)
  fthetaSetup = 0.;
        
  //Slits
  fSlitSize = 50.*mm;
  fSlit1Thickness = 5.*mm;
  fSlit2Thickness = 5.*mm;
  fSlit3Thickness = 5.*mm;
  fSlit4Thickness = 5.*mm;
  fSlit1SampleDistance = 350.*mm;
  fSlit2SampleDistance = 50.*mm;
  fSlit3SampleDistance = 50.*mm;
  fSlit4SampleDistance = 450.*mm;
  fSlit1xAperture = 4.*mm;
  fSlit2xAperture = 4.*mm;
  fSlit3xAperture = 4.*mm;
  fSlit4xAperture = 4.*mm;
  fSlit1yAperture = 4.*mm;
  fSlit2yAperture = 4.*mm;
  fSlit3yAperture = 4.*mm;
  fSlit4yAperture = 4.*mm;

  //Detector
  fDetectorSize = 20.*mm;
  fDetectorThickness = 20.*mm;
  fDetectorSampleDistance = 500.*mm;
        
  //Shielding
  fShieldingThickness = 4.*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SAXSDetectorConstruction::Construct()
{   
  //define custom materials
  DefineMaterials();
        
  //World    
  fWorldSolid = new G4Box("WorldSolid", 
                          fWorldSize*0.5, 
                          fWorldSize*0.5, 
                          fWorldSize*0.5);   
     
  fWorldLogic = new G4LogicalVolume(fWorldSolid, fAir, "WorldLogic");  
     
  fWorldLogic->SetVisAttributes(G4VisAttributes::GetInvisible());    
    
  fWorldPhysical = new G4PVPlacement(0, 
                                     G4ThreeVector(0., 0., 0.), 
                                     fWorldLogic, 
                                     "WorldPhysical", 
                                     0, 
                                     false, 
                                     0);
      
  //choose the phantom material
  switch (fPhantomMaterialIndex) {        
  case (1):                           
    fPhantomMaterial = fWater;
    break;
  case (2):                                           
    fPhantomMaterial = fMedMat;
    break;
  case (3):
    fPhantomMaterial = fPMMA;
    break;    
  case (4):
    fPhantomMaterial = fAdipose;
    break;   
  case (5):
    fPhantomMaterial = fGlandular;
    break; 
  case (6):
    fPhantomMaterial = fBreast5050;
    break; 
  case (7):
    fPhantomMaterial = fcarcinoma;
    break; 
  case (8):
    fPhantomMaterial = fkidney;
    break; 
  case (9):
    fPhantomMaterial = fliver;
    break; 
  case (10):
    fPhantomMaterial = fFat;
    break;  
  case (11):
    fPhantomMaterial = fBoneMatrix;
    break;  
  case (12):
    fPhantomMaterial = fMineral;
    break;  
  case (13):
    fPhantomMaterial = fBone;
    break;  
  case (14):
    fPhantomMaterial = ffatLowX;
    break; 
  case (15):
    fPhantomMaterial = fbonematrixLowX;
    break;
  case (16):
    fPhantomMaterial = fdryBoneLowX;
    break;
  case (17):
    fPhantomMaterial = fLexan;
    break;  
  case (18):
    fPhantomMaterial = fKapton;
    break;  
  case (19):
    fPhantomMaterial = fNylon;
    break;  
  case (20):
    fPhantomMaterial = fPolyethylene;
    break;  
  case (21):
    fPhantomMaterial = fPolystyrene;
    break; 
  case (22):
    fPhantomMaterial = fFormaline;
    break; 
  case (23):
    fPhantomMaterial = fAcetone;
    break; 
  case (24):
    fPhantomMaterial = fHperoxide;
    break;  
  case (25):
    fPhantomMaterial = fCIRS3070;
    break; 
  case (26):
    fPhantomMaterial = fCIRS5050;
    break; 
  case (27):
    fPhantomMaterial = fCIRS7030;
    break; 
  case (28):
    fPhantomMaterial = fRMI454;                          
    break;
  case (29):
    fPhantomMaterial = fAir;
    break; 
  case (30):
    fPhantomMaterial = fCustomMat;
    break;
  }    
        
  //Phantom (cylinder with axis orthogonal to the X-ray beam axis)
  G4Tubs* PhantomSolid = new G4Tubs("PhantomSolid", 
                                    0., 
                                    fPhantomDiameter*0.5, 
                                    fPhantomHeight*0.5, 
                                    0.*deg, 
                                    360.*deg);  
                                                        
  fPhantomLogic = new G4LogicalVolume(PhantomSolid, 
                                      fPhantomMaterial, 
                                      "PhantomLogic");        

  G4VisAttributes* PhantomVisAttribute =
    new G4VisAttributes(G4Colour(1., 1., 1.));  
  PhantomVisAttribute->SetForceSolid(true);
  fPhantomLogic->SetVisAttributes(PhantomVisAttribute);

  G4cout << "Phantom material: " << fPhantomMaterial->GetName() << G4endl;
  G4cout << "Phantom density: " << fPhantomMaterial->GetDensity()/(g/cm3) 
         << " g/cm3" << G4endl;
  G4cout << "Phantom mass: " << fPhantomLogic->GetMass()/g << " g" << G4endl;        

  G4double rotAngle = 90.*CLHEP::deg;
  G4RotationMatrix* PhantomRotationMatrix =
    new G4RotationMatrix(0., 0., 0.);
  PhantomRotationMatrix->rotateX(rotAngle);
        
  fPhantomPhysical = new G4PVPlacement(PhantomRotationMatrix, 
                                       G4ThreeVector(0., 0., fPhantomZ),
                                       fPhantomLogic,
                                       "PhantomPhysical",
                                       fWorldLogic,
                                       false,
                                       0);
          
  //setup rotation matrix (downstream of the phantom/sample)
  G4RotationMatrix* SetupRotationMatrix = new G4RotationMatrix();
  SetupRotationMatrix->rotateY(-fthetaSetup);
        
  //Slits
  G4Box* Slit1OutSolid = new G4Box("Slit1OutSolid",
                                   fSlitSize*0.5, 
                                   fSlitSize*0.5, 
                                   fSlit1Thickness*0.5); 
                                         
  G4Box* Slit2OutSolid = new G4Box("Slit2OutSolid", 
                                   fSlitSize*0.5, 
                                   fSlitSize*0.5,
                                   fSlit2Thickness*0.5);
                                         
  G4Box* Slit3OutSolid = new G4Box("Slit3OutSolid", 
                                   fSlitSize*0.5,
                                   fSlitSize*0.5, 
                                   fSlit3Thickness*0.5);
                                         
  G4Box* Slit4OutSolid = new G4Box("Slit4OutSolid", 
                                   fSlitSize*0.5, 
                                   fSlitSize*0.5, 
                                   fSlit4Thickness*0.5); 
         
  G4Box* Hole1Solid = new G4Box("Hole1Solid", 
                                fSlit1xAperture*0.5, 
                                fSlit1yAperture*0.5, 
                                fSlit1Thickness*0.51); 
                                      
  G4Box* Hole2Solid = new G4Box("Hole2Solid", 
                                fSlit2xAperture*0.5, 
                                fSlit2yAperture*0.5, 
                                fSlit2Thickness*0.51);  
                                      
  G4Box* Hole3Solid = new G4Box("Hole3Solid", 
                                fSlit3xAperture*0.5, 
                                fSlit3yAperture*0.5, 
                                fSlit3Thickness*0.51); 
                                       
  G4Box* Hole4Solid = new G4Box("Hole4Solid", 
                                fSlit4xAperture*0.5, 
                                fSlit4yAperture*0.5, 
                                fSlit4Thickness*0.51); 
        
  G4SubtractionSolid* Slit1Solid = new G4SubtractionSolid("Slit1Solid", 
                                                          Slit1OutSolid, 
                                                          Hole1Solid);
                                                                
  G4SubtractionSolid* Slit2Solid = new G4SubtractionSolid("Slit1Solid",
                                                          Slit2OutSolid,
                                                          Hole2Solid);
                                                                
  G4SubtractionSolid* Slit3Solid = new G4SubtractionSolid("Slit3Solid",
                                                          Slit3OutSolid,
                                                          Hole3Solid);
                                                                
  G4SubtractionSolid* Slit4Solid = new G4SubtractionSolid("Slit4Solid",
                                                          Slit4OutSolid, 
                                                          Hole4Solid);
     
  fSlit1Logic = new G4LogicalVolume(Slit1Solid, fTungsten, "Slit1Logic"); 
  fSlit2Logic = new G4LogicalVolume(Slit2Solid, fTungsten, "Slit2Logic");
  fSlit3Logic = new G4LogicalVolume(Slit3Solid, fTungsten, "Slit3Logic");
  fSlit4Logic = new G4LogicalVolume(Slit4Solid, fTungsten, "Slit4Logic");

  if (fIWantSlits) {
    G4cout << "Slit material: Tungsten" << G4endl;
    G4cout << "Slit1 thickness: " << fSlit1Thickness/mm << " mm" << G4endl;
    G4cout << "Slit2 thickness: " << fSlit2Thickness/mm << " mm" << G4endl;
    G4cout << "Slit3 thickness: " << fSlit3Thickness/mm << " mm" << G4endl;
    G4cout << "Slit4 thickness: " << fSlit4Thickness/mm << " mm" << G4endl;
    G4cout << "Slit1 aperture: " << fSlit1xAperture/mm << " x " 
           << fSlit1yAperture/mm << " mm2" << G4endl;
    G4cout << "Slit2 aperture: " << fSlit2xAperture/mm << " x " 
           << fSlit2yAperture/mm << " mm2" << G4endl;
    G4cout << "Slit3 aperture: " << fSlit3xAperture/mm << " x " 
           << fSlit3yAperture/mm << " mm2" << G4endl;
    G4cout << "Slit4 aperture: " << fSlit4xAperture/mm << " x " 
           << fSlit4yAperture/mm << " mm2" << G4endl;
  }
    
  G4VisAttributes* SlitlVisAttribute =
    new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  
  SlitlVisAttribute->SetForceSolid(true);
  fSlit1Logic->SetVisAttributes(SlitlVisAttribute);
  fSlit2Logic->SetVisAttributes(SlitlVisAttribute);
  fSlit3Logic->SetVisAttributes(SlitlVisAttribute);
  fSlit4Logic->SetVisAttributes(SlitlVisAttribute);
        
  G4double Slit1z = fPhantomZ - fSlit1SampleDistance;
  G4ThreeVector Slit1PositionVector = G4ThreeVector(0., 0., Slit1z); 
        
  G4double Slit2z = fPhantomZ - fSlit2SampleDistance;
  G4ThreeVector Slit2PositionVector = G4ThreeVector(0., 0., Slit2z);  
            
  G4double Slit3x = fSlit3SampleDistance*std::sin(fthetaSetup);
  G4double Slit3z = fPhantomZ + fSlit3SampleDistance*std::cos(fthetaSetup);
  G4ThreeVector Slit3PositionVector = G4ThreeVector(Slit3x, 0., Slit3z); 
    
  G4double Slit4x = fSlit4SampleDistance*std::sin(fthetaSetup);
  G4double Slit4z = fPhantomZ + fSlit4SampleDistance*std::cos(fthetaSetup);
  G4ThreeVector Slit4PositionVector = G4ThreeVector(Slit4x, 0., Slit4z);       
        
  if (fIWantSlits) {
    fSlit1Physical = new G4PVPlacement(0, 
                                       Slit1PositionVector,
                                       fSlit1Logic, 
                                       "Slit1Physical", 
                                       fWorldLogic, 
                                       false, 
                                       0);
                                          
    fSlit2Physical = new G4PVPlacement(0, 
                                       Slit2PositionVector, 
                                       fSlit2Logic, 
                                       "Slit2Physical", 
                                       fWorldLogic, 
                                       false, 
                                       0);
                                          
    fSlit3Physical = new G4PVPlacement(SetupRotationMatrix, 
                                       Slit3PositionVector, 
                                       fSlit3Logic, 
                                       "Slit3Physical", 
                                       fWorldLogic, 
                                       false, 
                                       0);
                                              
    fSlit4Physical = new G4PVPlacement(SetupRotationMatrix, 
                                       Slit4PositionVector, 
                                       fSlit4Logic, 
                                       "Slit4Physical", 
                                       fWorldLogic, 
                                       false, 
                                       0);
  }
                              
  //Detector (with shielding)
  G4Tubs* DetectorSolid = new G4Tubs("DetectorSolid", 
                                     0., 
                                     fDetectorSize*0.5, 
                                     fDetectorThickness*0.5, 
                                     0.*deg, 
                                     360.*deg); 
                
  fDetectorLogic = new G4LogicalVolume(DetectorSolid, fGe, "DetectorLogic");
        
  G4VisAttributes* DetectorVisAttribute =
    new G4VisAttributes(G4Colour(0., 0.5, 0.));  
  DetectorVisAttribute->SetForceSolid(true);
  fDetectorLogic->SetVisAttributes(DetectorVisAttribute);
    
  G4double Detx = fDetectorSampleDistance*std::sin(fthetaSetup);
  G4double Detz = fPhantomZ + fDetectorSampleDistance*std::cos(fthetaSetup);
  G4ThreeVector DetectorPositionVector = G4ThreeVector(Detx, 0., Detz); 
    
  fDetectorPhysical = new G4PVPlacement(SetupRotationMatrix, 
                                        DetectorPositionVector, 
                                        fDetectorLogic, 
                                        "DetectorPhysical", 
                                        fWorldLogic, 
                                        false, 
                                        0); 
           
  //Shielding
  G4double ShieldingSize = fDetectorSize+2*fShieldingThickness;
        
  G4double margin = 2.;
  G4double ShieldingLength = fDetectorThickness+fShieldingThickness*margin;        
        
  G4double ShieldingSampleDistance = fDetectorSampleDistance+
    fDetectorThickness*0.5+
    fShieldingThickness-
    ShieldingLength*0.5; 
                
  G4Tubs* ShieldingSolid = new G4Tubs("ShieldingSolid", 
                                      fDetectorSize*0.5, 
                                      ShieldingSize*0.5, 
                                      ShieldingLength*0.5, 
                                      0.*deg, 
                                      360.*deg); 
                
  fShieldingLogic = new G4LogicalVolume(ShieldingSolid, 
                                        fLead, 
                                        "ShieldingLogic");
        
  G4cout << "Shielding material: Lead" << G4endl;
  G4cout << "Shielding thickness: " << fShieldingThickness/mm 
         << " mm" << G4endl;
        
  G4VisAttributes* ShieldingVisAttribute =
    new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));  
  ShieldingVisAttribute->SetForceSolid(true);
  fShieldingLogic->SetVisAttributes(ShieldingVisAttribute);
        
  G4double Shieldx = ShieldingSampleDistance*std::sin(fthetaSetup);
  G4double Shieldz = fPhantomZ +
    ShieldingSampleDistance*std::cos(fthetaSetup);
  G4ThreeVector ShieldingPositionVector =
    G4ThreeVector(Shieldx, 0., Shieldz); 
   
  G4double ShieldingBackSampleDistance = fDetectorSampleDistance+
    fDetectorThickness*0.5+
    fShieldingThickness*0.5;
    
  G4Tubs* ShieldingBackSolid = new G4Tubs("ShieldingBackSolid", 
                                          0., 
                                          fDetectorSize*0.5, 
                                          fShieldingThickness*0.5, 
                                          0.*deg, 
                                          360.*deg); 
                                            
  fShieldingBackLogic = new G4LogicalVolume(ShieldingBackSolid, 
                                            fLead, 
                                            "ShieldingBackLogic");
    
  fShieldingBackLogic->SetVisAttributes(ShieldingVisAttribute);
   
  G4double ShieldBackx = ShieldingBackSampleDistance*std::sin(fthetaSetup);
  G4double ShieldBackz = fPhantomZ +
    ShieldingBackSampleDistance*std::cos(fthetaSetup);
  G4ThreeVector ShieldingBackPositionVector = 
    G4ThreeVector(ShieldBackx, 0., ShieldBackz);
         
  fShieldingPhysical = new G4PVPlacement(SetupRotationMatrix, 
                                         ShieldingPositionVector, 
                                         fShieldingLogic, 
                                         "ShieldingPhysical", 
                                         fWorldLogic, 
                                         false, 
                                         0);
                                              
  fShieldingBackPhysical = new G4PVPlacement(SetupRotationMatrix, 
                                             ShieldingBackPositionVector, 
                                             fShieldingBackLogic, 
                                             "ShieldingBackPhysical", 
                                             fWorldLogic, 
                                             false, 
                                             0);

  return fWorldPhysical;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSDetectorConstruction::ConstructSDandField()
{
  //Sensitive Volume
  G4VSensitiveDetector* vDetector = new SAXSSensitiveDetector("det");
  G4SDManager::GetSDMpointer()->AddNewDetector(vDetector); 
  fSensitiveVolume = fDetectorLogic;                               
  fSensitiveVolume->SetSensitiveDetector(vDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
