// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hTestStoppingPower.cc,v 1.3 2000-09-04 14:16:18 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
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
#include "G4hMollereNuclear.hh"


#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4Ions.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

HepTupleManager* hbookManager;

main()
{
  // ---- HBOOK initialization

  hbookManager = new HBookFile("stopping.paw", 68);
  assert (hbookManager != 0);
  
  // ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<G4endl;
    
  //--------- Materials definition ---------

  G4Material* Be = new G4Material("Beryllium", 4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6.,12.0*g/mole,2.265*g/cm3);
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* Si  = new G4Material("Silicon", 14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
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

  G4Material* water = new G4Material ("Water" ,"H_2O", 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" ,"C_2H_6", 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , "CsI", 4.53*g/cm3, 2);
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
  G4int numOfMaterials = theMaterialTable->length();

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

  G4double ecut = 10.0*mm;
  G4double pcut = 0.01*mm;
  electron->SetCuts(ecut);
  proton->SetCuts(pcut);
  antiproton->SetCuts(pcut);
  G4cout << "Cuts are following: cutElectron = " << ecut 
         << " mm; cutProton = " << pcut << " mm" << G4endl;  
  
  //--------- Ionisation processes definition and build physics table ------
    
  // Define models for parametrisation of electronic energy losses
  //  G4VLowEnergyModel* theBetheBlochModel = 
  //                   new G4hBetheBlochModel("Bethe-Bloch") ;
  //  theProtonModel = new G4hParametrisedLossModel(theProtonTable) ;
  //  theAntiProtonModel = new G4QAOLowEnergyLoss(theAntiProtonTable) ;
  //  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
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

  //  hIon[i]->SetNuclearStoppingOn();
  //  hIon[i]->SetStoppingPowerTableName("ICRU_R49p"); 
  
    theProcessManager[i]->AddProcess(hIon[i]);
    hIon[i]->BuildPhysicsTable(*part[i]);
  }


    
  //----------- Histogram -------------------------------------      

  G4cout << "Fill Hbook!" << G4endl;

  G4Material* material ;
 
  G4double minE = 1.0*eV, maxE = 10000.0*MeV, s;
  const G4int num = 200;
  G4double tkin;
  s = (log10(maxE)-log10(minE))/num;

 HepHistogram* h[60] ;
       
  // Test on Stopping Powers for all elements

 h[1] = hbookManager->histogram("p 40 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[2] =  hbookManager->histogram("p 100 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[3] =  hbookManager->histogram("p 400 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
    
 h[4] =  hbookManager->histogram("p 1 MeV   (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[5] =  hbookManager->histogram("p 4 MeV   (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;

 h[6] =  hbookManager->histogram("p 40 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[7] =  hbookManager->histogram("p 100 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[8] =  hbookManager->histogram("p 400 keV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[9] =  hbookManager->histogram("p 1 MeV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;
 h[10] =  hbookManager->histogram("p 4 MeV   (keV*cm2/10^15!atoms) ICRU_49p"
                                                  ,92,0.5,92.5) ;

 h[11] =  hbookManager->histogram("p 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[12] =  hbookManager->histogram("p 100 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[13] =  hbookManager->histogram("p 400 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[14] =  hbookManager->histogram("p 1 MeV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[15] =  hbookManager->histogram("p 4 MeV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;

 h[16] =  hbookManager->histogram("He 10 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[17] =  hbookManager->histogram("He 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[18] =  hbookManager->histogram("He 100 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 
 h[19] =  hbookManager->histogram("He 400 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[20] =  hbookManager->histogram("He 1 MeV  (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[21] =  hbookManager->histogram("He 4 MeV   (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;
 h[22] =  hbookManager->histogram("He 10 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[23] =  hbookManager->histogram("He 40 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[24] =  hbookManager->histogram("He 100 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[25] =  hbookManager->histogram("He 400 keV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[26] =  hbookManager->histogram("He 1 MeV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;
 h[27] =  hbookManager->histogram("He 4 MeV   (keV*cm2/10^15!atoms) ICRU49He"
                                                  ,92,0.5,92.5) ;

 h[28]= hbookManager->histogram("p   in C (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[29]= hbookManager->histogram("p   in Al (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[30]= hbookManager->histogram("p   in Si (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[31]= hbookManager->histogram("p   in Cu (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[32]= hbookManager->histogram("p   in Fe (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[33]= hbookManager->histogram("p   in Pb (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[34]= hbookManager->histogram("p   in C2H6 (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[35]= hbookManager->histogram("p   in H2O (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[36]= hbookManager->histogram("p   in lAr (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[37]= hbookManager->histogram("p   in CsI (MeV/mm) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;


 h[38]= hbookManager->histogram("p   in C (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[39]= hbookManager->histogram("p   in Al (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[40]= hbookManager->histogram("p   in Si (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[41]= hbookManager->histogram("p   in Cu (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[42]= hbookManager->histogram("p   in Fe (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[43]= hbookManager->histogram("p   in Pb (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[44]= hbookManager->histogram("p   in C2H6 (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[45]= hbookManager->histogram("p   in H2O (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[46]= hbookManager->histogram("p   in lAr (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[47]= hbookManager->histogram("p   in CsI (MeV/mm)Ziegler1977p"
                                   ,num,log10(minE),log10(maxE)) ;


 h[48]= hbookManager->histogram("He effective charge for Cu"
                                   ,num,log10(minE),log10(maxE)) ;
 h[49]= hbookManager->histogram("C12 effective charge in Cu"
                                   ,num,log10(minE),log10(maxE)) ;

 h[50]= hbookManager->histogram("He in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[51]= hbookManager->histogram("C12 in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;
 h[52]= hbookManager->histogram("Ar40 in Al (MeV/(mg/cm2)) ICRU49p"
                                   ,num,log10(minE),log10(maxE)) ;

 // Table with the data
 h[53] = hbookManager->histogram("Data p 40 keV (keV*cm2/10^15!atoms) Ziegler1977p"
                                                  ,92,0.5,92.5) ;
 h[54] =  hbookManager->histogram("Data He 40 keV (keV*cm2/10^15!atoms) Ziegler1977He"
                                                  ,92,0.5,92.5) ;

 h[55] = hbookManager->histogram("p 40 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[56] =  hbookManager->histogram("p 100 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[57] =  hbookManager->histogram("p 400 keV (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;    
 h[58] =  hbookManager->histogram("p 1 MeV   (keV*cm2/10^15!atoms) Ziegler1985p"
                                                  ,92,0.5,92.5) ;
 h[59] =  hbookManager->histogram("p 4 MeV   (keV*cm2/10^15!atoms) Ziegler1985p"
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
    h[1]->accumulate(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[55]->accumulate(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[6]->accumulate(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[11]->accumulate(z,de) ;
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[17]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[23]->accumulate(z,de) ;
  }
  tau = 100.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[2]->accumulate(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[56]->accumulate(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[7]->accumulate(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[12]->accumulate(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[18]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[24]->accumulate(z,de) ;
  }
  
  tau = 400.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[3]->accumulate(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[57]->accumulate(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[8]->accumulate(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[13]->accumulate(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[19]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[25]->accumulate(z,de) ;
  }
  tau = 1000.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[4]->accumulate(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[58]->accumulate(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[9]->accumulate(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[14]->accumulate(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[20]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[26]->accumulate(z,de) ;
  }
  tau = 4000*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    de = Z77p->ElectronicStoppingPower(z, tau) ;
    h[5]->accumulate(z,de) ;
    de = Z85p->ElectronicStoppingPower(z, tau) ;
    h[59]->accumulate(z,de) ;
    de = I49p->ElectronicStoppingPower(z, tau) ;
    h[10]->accumulate(z,de) ;
    de = Z77He->ElectronicStoppingPower(z, tau) ;
    h[15]->accumulate(z,de) ;
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[21]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[27]->accumulate(z,de) ;
  }
  tau = 10.0*keV ;
  for ( j=1; j<93; j++)
  { 
    z = G4double(j);
    q2    = theIonEffChargeModel->TheValue(part[3],el[j],tkin);
    de = (Z77He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[16]->accumulate(z,de) ;
    de = (I49He->ElectronicStoppingPower(z, tau/rateMass))*q2 ;
    h[22]->accumulate(z,de) ;
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
    h[53]->accumulate(double(j),p40[j-1]) ;

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
    h[54]->accumulate(double(j),he40[j-1]) ;
  }

  G4cout << "Stopping Power histograms are filled!" << G4endl;

  // dedx
  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = hIon[0]->ComputeDEDX(part[0],Graphite,tkin) ;
    h[28]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Al,tkin) ;
    h[29]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Si,tkin) ;
    h[30]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Cu,tkin) ;
    h[31]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Fe,tkin) ;
    h[32]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],Pb,tkin) ;
    h[33]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],ethane,tkin) ;
    //    G4cout << "ethane: E = " << tkin << "; dedx = " << de << G4endl;
    h[34]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],water,tkin) ;
    //    G4cout << "water : E = " << tkin << "; dedx = " << de << G4endl;
    h[35]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],LAr,tkin) ;
    h[36]->accumulate(log10(tkin),de) ;
    de = hIon[0]->ComputeDEDX(part[0],csi,tkin) ;
    h[37]->accumulate(log10(tkin),de) ;
  }

  G4cout << "Proton's dEdx histograms are filled!" << G4endl;

  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = theIonEffChargeModel->TheValue(part[3],Cu,tkin) ;
    //  G4cout << "E = " << tkin << "; dedx = " << de << G4endl;
    h[48]->accumulate(log10(tkin),de) ;
    de = theIonEffChargeModel->TheValue(part[4],Cu,tkin) ;
    h[49]->accumulate(log10(tkin),de) ;
    de = hIon[3]->ComputeDEDX(part[3],Al,tkin) ;
    h[50]->accumulate(log10(tkin),de) ;
    de = hIon[4]->ComputeDEDX(part[4],Al,tkin) ;
    h[51]->accumulate(log10(tkin),de) ;
    de = hIon[5]->ComputeDEDX(part[5],Al,tkin) ;
    h[52]->accumulate(log10(tkin),de) ;
  }

  G4cout << "Ions's dEdx histograms are filled!" << G4endl;

  theProcessManager[7] = new G4ProcessManager(part[0]);
  part[0]->SetProcessManager(theProcessManager[7]);
  hIon[7] = new G4hLowEnergyIonisation();
  hIon[7]->SetElectronicStoppingPowerModel(part[0],"Ziegler1977p"); 
  theProcessManager[7]->AddProcess(hIon[7]);
  hIon[7]->BuildPhysicsTable(*part[0]);

  G4cout << "Ziegler's dEdx histograms will be filled" << G4endl;

  for (j = 0 ; j < num-1 ; j++) {
    tkin = pow(10.0,(log10(minE) + (G4double(j)+0.5)*s));
    de = hIon[7]->ComputeDEDX(part[0],Graphite,tkin) ;
    h[38]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Al,tkin) ;
    h[39]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Si,tkin) ;
    h[40]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Cu,tkin) ;
    h[41]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Fe,tkin) ;
    h[42]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],Pb,tkin) ;
    h[43]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],ethane,tkin) ;
    //    G4cout << "ethane: E = " << tkin << "; dedx = " << de << G4endl;
    h[44]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],water,tkin) ;
    //    G4cout << "water:  E = " << tkin << "; dedx = " << de << G4endl;
    h[45]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],LAr,tkin) ;
    h[46]->accumulate(log10(tkin),de) ;
    de = hIon[7]->ComputeDEDX(part[0],csi,tkin) ;
    h[47]->accumulate(log10(tkin),de) ;
  }
 
  //---------------------- Fill Ntuple ------------------------
 
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
				   
   
  //----------- End of work -------------------------------------      

  G4cout << "Save Ntuple and Hbook" << G4endl;  
  hbookManager->write();

  delete hbookManager;
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
