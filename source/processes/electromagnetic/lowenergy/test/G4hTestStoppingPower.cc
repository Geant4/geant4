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
// $Id: G4hTestStoppingPower.cc,v 1.13 2002-08-08 16:58:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4hTestStoppingPower.cc
//
//      Author:        Vladimir Ivanchenko
// 
//      Creation date: 28 July 2000
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ProcessManager.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hBetheBlochModel.hh"
#include "G4hParametrisedLossModel.hh"
#include "G4hNuclearStoppingModel.hh"
#include "G4QAOLowEnergyLoss.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4hZiegler1977p.hh"
#include "G4hZiegler1985p.hh"
#include "G4hZiegler1977He.hh"
#include "G4hICRU49p.hh"
#include "G4hICRU49He.hh"

#include "G4hZiegler1977Nuclear.hh"
#include "G4hZiegler1985Nuclear.hh"
//#include "G4hMollereNuclear.hh" // exist no more


#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4Ions.hh"

// New Histogramming (from AIDA and Anaphe):
#include <memory> // for the auto_ptr(T>

#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"

#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
//#include "AIDA/IHistogram3D.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "hTest/include/G4IonC12.hh"
#include "hTest/include/G4IonAr40.hh"

#include "G4Timer.hh"

main()
{
  // ---- HBOOK initialization

  G4String hFile = "htest.hbook";
    
  //--------- Materials definition ---------

  G4Material* Be = new G4Material("Beryllium", 4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6.,12.0*g/mole,2.265*g/cm3);
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* Si  = new G4Material("Silicon", 14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
  G4Material* Ti  = new G4Material("Titan", 22., 47.867*g/mole, 4.54*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten",74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",   82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  //  G4Material* water = new G4Material ("Water" ,"H_2O", 1.*g/cm3, 2);
  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  //  G4Material* ethane = new G4Material ("Ethane" ,"C_2H_6", 0.4241*g/cm3, 2);
  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  //  G4Material* csi = new G4Material ("CsI" , "CsI", 4.53*g/cm3, 2);
  G4Material* csi = new G4Material ("CsI" ,  4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  G4Material* el[93];
  G4double z;
  G4int j;

  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    el[j] = new G4Material ("atom",  z , 2.0*z*g/mole, 10.0 *g/cm3);
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  
  //  create table
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  G4double dimx = 1*mm, dimy = 1*mm, dimz = 1*mm;
  G4int imat = 0;

  G4cout<<"The material is: "<<(*theMaterialTable)(imat)->GetName()<<endl;

  // Geometry definitions
  G4Box* theFrame = new G4Box ("Frame",dimx, dimy, dimz);
  
  G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
				      (*theMaterialTable)(imat),
						      "LFrame", 0, 0, 0);
  
  G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
				       "PFrame",LogicalFrame,0,false,0);


  //--------- Particle definition ---------
  G4Electron* theElectron = G4Electron::Electron();
  
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiproton =G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* deuteron = G4Deuteron::DeuteronDefinition();
  G4ParticleDefinition* alpha = G4Alpha::AlphaDefinition();

  G4Ions* iC12 = new G4Ions::G4Ions(
              "IonC12",    11.14945*GeV,       0.0*MeV,  +6.0*eplus, 
		    0,              +1,             0,          
		    0,               0,             0,             
	    "static_nucleus",        0,            +12,           0,
		 true,            -1.0,          NULL);

  G4Ions* iAr40 = new G4Ions::G4Ions(
              "IonAr40",    37.291*GeV,       0.0*MeV,  +18.0*eplus, 
		    0,              +1,             0,          
		    0,               0,             0,             
	    "static_nucleus",        0,            +40,           0,
		 true,            -1.0,          NULL);

  G4Ions* iFe56 = new G4Ions::G4Ions(
              "IonFe56",    52.0308*GeV,       0.0*MeV,  +26.0*eplus, 
		    0,              +1,             0,          
		    0,               0,             0,             
	    "static_nucleus",        0,            +56,           0,
		 true,            -1.0,          NULL);

  G4ParticleDefinition* ionC12  = iC12->IonsDefinition();
  G4ParticleDefinition* ionAr40 = iAr40->IonsDefinition();
  G4ParticleDefinition* ionFe56 = iFe56->IonsDefinition();

  G4ParticleDefinition* part[7];
  part[0] = proton;  
  part[1] = antiproton;  
  part[2] = deuteron;  
  part[3] = alpha;  
  part[4] = ionC12;  
  part[5] = ionAr40;  
  part[6] = ionFe56;  

  G4double ecut = 1000.0*mm;
  G4double pcut = 0.0001*mm;
  electron->SetCuts(ecut);
  proton->SetCuts(pcut);
  antiproton->SetCuts(pcut);
  G4cout << "Cuts are following: cutElectron = " << ecut 
         << " mm; cutProton = " << pcut << " mm" << G4endl;  
  
  //--------- Ionisation processes definition and build physics table --
    
  // Define models for parametrisation of electronic energy losses
  //  G4VLowEnergyModel* theBetheBlochModel = 
  //                   new G4hBetheBlochModel("Bethe-Bloch") ;
  //  theProtonModel = new G4hParametrisedLossModel(theProtonTable) ;
  //  theAntiProtonModel = new G4QAOLowEnergyLoss(theAntiProtonTable) ;
  
  G4hNuclearStoppingModel* theNuclearStoppingModel = 
    //                    new G4hNuclearStoppingModel("Ziegler1977") ;
  //                     new G4hNuclearStoppingModel("Ziegler1985") ;
                       new G4hNuclearStoppingModel("ICRU_R49") ;

  G4VLowEnergyModel* theIonEffChargeModel = 
                     new G4hIonEffChargeSquare("Ziegler1988") ;

  
  G4hLowEnergyIonisation* hIon[8];
  G4ProcessManager* theProcessManager[8];
  G4int i;

  G4cout << "Define processes!" << G4endl;

  for( i=0; i<7; i++) {    
    G4cout << "Ionisation process for particle " << i << G4endl;
    part[i]->SetCuts(pcut);
    theProcessManager[i] = new G4ProcessManager(part[i]);
    part[i]->SetProcessManager(theProcessManager[i]);
    hIon[i] = new G4hLowEnergyIonisation();
    hIon[i]->SetEnlossFluc(false) ;

    // hIon[i]->SetBarkasOff();
    //    hIon[i]->SetNuclearStoppingOff();
  //  hIon[i]->SetStoppingPowerTableName("ICRU_R49p"); 
  
    theProcessManager[i]->AddProcess(hIon[i]);
    hIon[i]->BuildPhysicsTable(*part[i]);
  }


    
  //----------- Histogram -------------------------------------      

  G4cout << "Fill Hbook!" << G4endl;

    // Creating the analysis factory
    G4std::auto_ptr< IAnalysisFactory > af( AIDA_createAnalysisFactory() );

    // Creating the tree factory
    G4std::auto_ptr< ITreeFactory > tf( af->createTreeFactory() );

    // Creating a tree mapped to a new hbook file.
    G4std::auto_ptr< ITree > tree( tf->create( hFile,false,false,"hbook" ) );
    G4std::cout << "Tree store : " << tree->storeName() << G4std::endl;

  G4Material* material ;
 
  G4double minE = 1.0*eV, maxE = 10000.0*MeV, s;
  const G4int num = 200;
  G4double tkin;
  s = (log10(maxE)-log10(minE))/num;

  IHistogram1D* h[71] ;
       
  // Test on Stopping Powers for all elements

 h[1] = hf->create1D("1","p 40 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[2] =  hf->create1D("2","p 100 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[3] =  hf->create1D("3","p 400 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
    
 h[4] =  hf->create1D("4","p 1 MeV   (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[5] =  hf->create1D("5","p 4 MeV   (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;

 h[6] =  hf->create1D("6","p 40 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[7] =  hf->create1D("7","p 100 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[8] =  hf->create1D("8","p 400 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[9] =  hf->create1D("9","p 1 MeV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[10] =  hf->create1D("10","p 4 MeV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;

 h[11] =  hf->create1D("11","p 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[12] =  hf->create1D("12","p 100 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[13] =  hf->create1D("13","p 400 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[14] =  hf->create1D("14","p 1 MeV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[15] =  hf->create1D("15","p 4 MeV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;

 h[16] =  hf->create1D("16","He 10 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[17] =  hf->create1D("17","He 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[18] =  hf->create1D("18","He 100 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 
 h[19] =  hf->create1D("19","He 400 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[20] =  hf->create1D("20","He 1 MeV  (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[21] =  hf->create1D("21","He 4 MeV   (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[22] =  hf->create1D("22","He 10 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[23] =  hf->create1D("23","He 40 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[24] =  hf->create1D("24","He 100 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[25] =  hf->create1D("25","He 400 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[26] =  hf->create1D("26","He 1 MeV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[27] =  hf->create1D("27","He 4 MeV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;

 h[28]= hf->create1D("28","p   in C (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[29]= hf->create1D("29","p   in Al (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[30]= hf->create1D("30","p   in Si (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[31]= hf->create1D("31","p   in Cu (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[32]= hf->create1D("32","p   in Fe (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[33]= hf->create1D("33","p   in Pb (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[34]= hf->create1D("34","p   in C2H6 (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[35]= hf->create1D("35","p   in H2O (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[36]= hf->create1D("36","p   in lAr (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[37]= hf->create1D("37","p   in CsI (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;


 h[38]= hf->create1D("38","p   in C (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[39]= hf->create1D("39","p   in Al (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[40]= hf->create1D("40","p   in Si (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[41]= hf->create1D("41","p   in Cu (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[42]= hf->create1D("42","p   in Fe (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[43]= hf->create1D("43","p   in Pb (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[44]= hf->create1D("44","p   in C2H6 (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[45]= hf->create1D("45","p   in H2O (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[46]= hf->create1D("46","p   in lAr (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[47]= hf->create1D("47","p   in CsI (MeV/mm)Ziegler1985p"
                                   ,num,log10(minE),log10(maxE)) ;


 h[48]= hf->create1D("48","He effective charge for Cu"
                                   ,num,log10(minE),log10(maxE)) ;
 h[49]= hf->create1D("49","C12 effective charge in Cu"
                                   ,num,log10(minE),log10(maxE)) ;

 h[50]= hf->create1D("50","He in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[51]= hf->create1D("51","C12 in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[52]= hf->create1D("52","Ar40 in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;

 // Table with the data
 h[53] = hf->create1D("53","Data p 40 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[54] =  hf->create1D("54","Data He 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;

 h[55] = hf->create1D("55","p 40 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[56] =  hf->create1D("56","p 100 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[57] =  hf->create1D("57","p 400 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;    
 h[58] =  hf->create1D("58","p 1 MeV   (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[59] =  hf->create1D("59","p 4 MeV   (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
// Histo for Antiproton
 
 h[60]= hf->create1D("60","pbar   in C (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[61]= hf->create1D("61","pbar   in Al (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[62]= hf->create1D("62","pbar   in Si (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[63]= hf->create1D("63","pbar   in Cu (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[64]= hf->create1D("64","pbar   in Fe (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[65]= hf->create1D("65","pbar   in Pb (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[66]= hf->create1D("66","pbar   in C2H6 (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[67]= hf->create1D("67","pbar   in H2O (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[68]= hf->create1D("68","pbar   in lAr (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;
 h[69]= hf->create1D("69","pbar   in CsI (MeV/mm) QAOLoss"
                                   ,num,log10(minE),log10(maxE)) ;

 h[70] = hf->create1D("70","p 6.5 MeV (keV*cm2/10^15!atoms) Ziegler77p"
                                                  ,92,0.5,92.5) ;

 
 G4VhElectronicStoppingPower* Z77p = new G4hZiegler1977p() ;
 G4VhElectronicStoppingPower* Z85p = new G4hZiegler1985p() ;
 G4VhElectronicStoppingPower* Z77He = new G4hZiegler1977He() ;
 G4VhElectronicStoppingPower* I49p = new G4hICRU49p() ;
 G4VhElectronicStoppingPower* I49He = new G4hICRU49He() ;

  G4double de;

  // Ziegler1977p
  G4double tau = 40.0*keV ;
  G4double q2;
  G4double rateMass=4.0026/1.007276;

  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[1]->fill(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[55]->fill(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[6]->fill(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[11]->fill(z,de) ;
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[17]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[23]->fill(z,de) ;
  }
  tau = 100.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[2]->fill(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[56]->fill(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[7]->fill(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[12]->fill(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[18]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[24]->fill(z,de) ;
  }
  
  tau = 400.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[3]->fill(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[57]->fill(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[8]->fill(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[13]->fill(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[19]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[25]->fill(z,de) ;
  }
  tau = 1000.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[4]->fill(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[58]->fill(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[9]->fill(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[14]->fill(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[20]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[26]->fill(z,de) ;
  }
  tau = 4000*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[5]->fill(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[59]->fill(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[10]->fill(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[15]->fill(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[21]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[27]->fill(z,de) ;
  }
  tau = 6.5*MeV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[70]->fill(z,de) ;
  }
  tau = 10.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[16]->fill(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[22]->fill(z,de) ;
  }

  for ( j=1; j<93; j++)
  { 
    static G4double p40[92] = {
    6.22, 6.55, 7.61, 10.4, 12.9, 13.9, 15.9, 14.6, 11.7, 11.1,
    14.2, 20.6, 20.7,  22.3, 18.1, 19.3, 27.5, 30.5, 27.9, 29.8,
    28.3, 26.7, 25.1,  22.5, 19.8, 20.1, 18.0, 20.1, 20.4, 23.4,
    27.5, 29.0, 29.3,  31.0, 31.1, 34.8, 31.5, 35.0, 35.5, 37.3,
    38.3, 35.9, 38.0,  34.4, 33.4, 29.8, 30.9, 32.7, 34.9, 36.0,
    40.9, 39.1, 43.1,  45.8, 40.8, 44.1, 44.9, 42.0, 41.0, 40.0,
    39.0, 38.0, 37.1,  38.2, 35.3, 31.5, 29.9, 29.1, 28.3, 27.5, 
    28.1, 28.9, 27.3,  26.3, 29.9, 29.1, 28.4, 25.8, 28.0, 24.9,
    27.3, 30.7, 34.3,  35.4, 35.6, 35.5, 39.8, 42.9, 43.7, 44.1,
    42.4, 41.8} ;
    h[53]->fill(double(j),p40[j-1]) ;

    static G4double he40[92] = {
    11.8, 16.7, 21.9,  24.1, 34.7, 36.1, 45.5, 46.1, 45.8, 44.9,
    46.9, 51.7, 54.9,  63.2, 63.8, 67.4, 79.8, 80.8, 78.0, 83.2,
    85.4, 84.2, 90.6,  85.7, 84.1, 86.4, 82.0, 77.7, 74.4, 75.1,
    75.6, 80.7, 84.9,  87.0, 88.6, 103.0, 103.0, 113.0, 116.0, 123.0,
    120.0, 115.0, 118.0, 115.0, 113.0, 112.0, 111.0, 110.0, 116.0, 117.0,
    119.0, 123.0, 122.0, 139.0, 138.0, 143.0, 152.0, 141.0, 139.0, 137.0,
    135.0, 134.0, 130.0, 137.0, 129.0, 128.0, 124.0, 125.0, 120.0, 118.0, 
    119.0, 120.0, 121.0, 122.0, 125.0, 127.0, 128.0, 120.0, 125.0, 128.0,
    136.0, 138.0, 144.0, 148.0, 151.0, 162.0, 166.0, 177.0, 179.0, 183.0, 
    177.0, 176.0 } ;   
    h[54]->fill(double(j),he40[j-1]) ;
  }

  G4cout << "Stopping Power histograms are filled!" << G4endl;

  // dedx
  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = hIon[0]->ComputeDEDX(part[0],Graphite,tkin) ;
    h[28]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Al,tkin) ;
    h[29]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Si,tkin) ;
    h[30]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Cu,tkin) ;
    h[31]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Fe,tkin) ;
    h[32]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Pb,tkin) ;
    h[33]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],ethane,tkin) ;
    //    G4cout << "ethane: E = " << tkin << "; dedx = " << de << G4endl;
    h[34]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],water,tkin) ;
    //    de += theNuclearStoppingModel->TheValue(part[1],water,tkin) ;
    //    G4cout << "water : E = " << tkin << "; dedx = " << de << G4endl;
    h[35]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],LAr,tkin) ;
    h[36]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],csi,tkin) ;
    h[37]->fill(log10(tkin),de) ;
  }

  G4cout << "Proton's dEdx histograms are filled!" << G4endl;
  
  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = hIon[0]->ComputeDEDX(part[1],Graphite,tkin) ;
    h[60]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],Al,tkin) ;
    h[61]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],Si,tkin) ;
    h[62]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],Cu,tkin) ;
    h[63]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],Fe,tkin) ;
    h[64]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],Pb,tkin) ;
    h[65]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],ethane,tkin) ;
    //    G4cout << "ethane: E = " << tkin << "; dedx = " << de << G4endl;
    h[66]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],water,tkin) ;
    // G4cout << "water : E = " << tkin << "; dedx = " << de << G4endl;
    h[67]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],LAr,tkin) ;
    h[68]->fill(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[1],csi,tkin) ;
    h[69]->fill(log10(tkin),de) ;
  }

  G4cout << "AntiProton's dEdx histograms are filled!" << G4endl;

  G4double mProt = part[0]->GetPDGMass()*1.007276;
  G4double fact  = cm/(2700.0*MeV) ;      // to MeV/mg/cm^2

  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = theIonEffChargeModel->TheValue(part[3],Cu,tkin) ;
    //  G4cout << "E = " << tkin << "; dedx = " << de << G4endl;
    h[48]->fill(log10(tkin),de) ;
    de = theIonEffChargeModel->TheValue(part[4],Cu,tkin) ;
    h[49]->fill(log10(tkin),de) ;
    G4double tRed = tkin * (part[3]->GetPDGMass())/mProt ;
    de = hIon[3]->ComputeDEDX(part[3],Al,tRed) ;
    de += theNuclearStoppingModel->TheValue(part[3],Al,tRed) ;
    h[50]->fill(log10(tkin),de*fact) ;
    tRed = tkin * (part[4]->GetPDGMass())/mProt ;
    de = hIon[4]->ComputeDEDX(part[4],Al,tRed) ;
    //de += theNuclearStoppingModel->TheValue(part[4],Al,tRed) ;
    h[51]->fill(log10(tkin),de*fact) ;
    tRed = tkin * (part[5]->GetPDGMass())/mProt ;
    de = hIon[5]->ComputeDEDX(part[5],Al,tRed) ;
    //de += theNuclearStoppingModel->TheValue(part[5],Al,tRed) ;
    h[52]->fill(log10(tkin),de*fact) ;
  } 

  G4cout << "Ions's dEdx histograms are filled!" << G4endl;

  theProcessManager[7] = new G4ProcessManager(part[0]);
  part[0]->SetProcessManager(theProcessManager[7]);
  hIon[7] = new G4hLowEnergyIonisation();
  hIon[7]->SetElectronicStoppingPowerModel(part[0],"Ziegler1985p"); 
  hIon[7]->SetEnlossFluc(false) ;
  //hIon[7]->SetNuclearStoppingOff();
  //  hIon[7]->SetBarkasOff();
  theProcessManager[7]->AddProcess(hIon[7]);
  hIon[7]->BuildPhysicsTable(*part[0]);

  G4cout << "Ziegler's dEdx histograms will be filled" << G4endl;

  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = hIon[7]->ComputeDEDX(part[0],Graphite,tkin) ;
    h[38]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Al,tkin) ;
    h[39]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Si,tkin) ;
    h[40]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Cu,tkin) ;
    h[41]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Fe,tkin) ;
    h[42]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Pb,tkin) ;
    h[43]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],ethane,tkin) ;
    //    G4cout << "ethane: E = " << tkin << "; dedx = " << de << G4endl;
    h[44]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],water,tkin) ;
    //    G4cout << "water:  E = " << tkin << "; dedx = " << de << G4endl;
    h[45]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],LAr,tkin) ;
    h[46]->fill(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],csi,tkin) ;
    h[47]->fill(log10(tkin),de) ;
  }
 
  //---------------------- Fill Ntuple ------------------------
  /* 
  G4cout << "Fill Ntuple!" << G4endl;

  // ---- primary ntuple ------
  HepTuple* ntuple = hbookManager->ntuple("dEdx ntuple");
  assert (ntuple != 0);
  
  for( i = 0; i < num; i++) { 
    tkin = pow(10,(log10(minE) + i*s));
    
    for ( G4int j = 0 ; j < numOfMaterials; j++ ) {

      material = (*theMaterialTable)[j] ;
      // get elements in the actual material,
      const G4ElementVector* theElementVector = material->GetElementVector() ;
      const G4double* theAtomicNumDensityVector = 
                         material->GetAtomicNumDensityVector() ;
      const G4int NumberOfElements = material->GetNumberOfElements() ;
  
      //  loop for the elements in the material
      //  to find out average values Z, vF, lF
      G4double z = 0.0, norm = 0.0 ; 

      if( 1 == NumberOfElements ) {
        z = material->GetZ() ;

      } else {
        for (G4int iel=0; iel<NumberOfElements; iel++) {
          const G4Element* element = (*theElementVector)(iel) ;
          G4double z2 = element->GetZ() ;
          const G4double weight = theAtomicNumDensityVector[iel] ;
          norm += weight ;
          z    += z2 * weight ;
        }
        z  /= norm ;
      }
   
      for (G4int k=0 ; k<7; k++) {
         G4double dedx = hIon[k]->ComputeDEDX(part[k],material,tkin) ;
         G4double q = theIonEffChargeModel->TheValue(part[k],material,tkin);
         ntuple->column("part",k);
         ntuple->column("mat",j);
         ntuple->column("ie",i);
         ntuple->column("zeff",z);
         ntuple->column("tkin",tkin/MeV);
         ntuple->column("mass",(part[k]->GetPDGMass())/MeV);
         ntuple->column("char",(part[k]->GetPDGCharge())/eplus);
         ntuple->column("cha2",q);
         ntuple->column("dedx",dedx*mm/MeV);
         ntuple->dumpData();
      }
    }
  }
  */				      
  //----------- End of work -------------------------------------      

      G4std::cout << "Committing..." << G4std::endl;
      tree->commit();
      G4std::cout << "Closing the tree..." << G4std::endl;
      tree->close();

  G4cout << "Ntuple and Hbook are saved" << G4endl;

  // delete materials and elements
  delete Be;
  delete Graphite;
  delete Al;
  delete Si;
  delete LAr;
  delete Fe;
  delete Cu;
  delete W;
  delete Pb;
  delete U;
  delete H;
  delete C;
  delete Cs;
  delete I;
  delete O;
  delete water;
  delete ethane;
  delete csi;
  G4cout << "Materials are deleted" << G4endl;
  delete theIonEffChargeModel;
  delete Z77p, Z77He, I49p, I49He;

  cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}  
