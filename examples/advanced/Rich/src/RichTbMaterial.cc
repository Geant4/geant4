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
// Rich advanced example for Geant4
// RichTbMaterial.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include "globals.hh"
#include "RichTbMaterial.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "RichTbMaterialParameters.hh"
#include "RichTbGeometryParameters.hh"
#include "G4MaterialPropertyVector.hh"


RichTbMaterial::RichTbMaterial(RichTbRunConfig* RConfig):
   RichTbAerogelMaterial(std::vector<G4Material*> (MaxNumberOfAerogelTypes)),
   RichTbFilterMaterial(std::vector<G4Material*>(MaxNumberOfFilterTypes)){ 

  rConfig=RConfig;

  G4double a,z,density;  //a=mass of a mole;
  // z=mean number of protons;
  G4String name,symbol;
  //G4int isz,isn;    //isz= number of protons in an isotope;
  //isn= number of nucleons in an isotope;
  
  
  G4int numel,natoms;  //numel=Number of elements constituting a material.                     
  G4double fractionmass;
  G4double temperature, pressure;
  // G4double FactorOne=1.0;
  G4UnitDefinition::BuildUnitsTable();
  
 //PhotonEnergy
  G4int ibin=0;
  G4double PhotonEnergyStep=(PhotonMaxEnergy-PhotonMinEnergy)/
                            NumPhotWaveLengthBins;
  G4double* PhotonMomentum=new G4double[NumPhotWaveLengthBins];
  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    PhotonMomentum[ibin]=PhotonMinEnergy+PhotonEnergyStep*ibin;
  }

  G4cout << "\nNow Define Elements ..\n" <<G4endl;

// Nitrogen

 a=14.01*g/mole;
 G4Element* elN = new G4Element(name="Nitrogen", 
                               symbol="N", z=7., a);

//Oxygen

 a=16.00*g/mole;
 G4Element* elO = new G4Element(name="Oxygen", 
                               symbol="O", z=8., a);

//Hydrogen

 a=1.01*g/mole;
 G4Element* elH = new G4Element(name="Hydrogen",
                               symbol="H",z=1.,a);

//Carbon

 a=12.01*g/mole;
 G4Element* elC = new G4Element(name="Carbon", 
                               symbol="C",z=6.,a);

//Silicon

 a=28.09*g/mole;
 G4Element* elSi = new G4Element(name="Silicon",
                                symbol="Si",z=14.,a);  
 //Fluorine
 a=18.998*g/mole;
 G4Element* elF = new G4Element(name="Fluorine",
                                symbol="F",z=9.,a);  
 //Aluminum
 a=26.98*g/mole;
 G4Element* elAL =new G4Element(name="Aluminium",
                                symbol="Al",z=13.,a); 

 //Sodium
 a=22.99*g/mole; 
 G4Element* elNa = new G4Element(name="Sodium",
                                symbol="Na",z=11.,a);  

 //Potassium
 a=39.10*g/mole; 
 G4Element* elK = new G4Element(name="Potassium",
                                symbol="K",z=19.,a);  

 //Cesium

 // a=132.91*g/mole; 
 // G4Element* elCs = new G4Element(name="Cesium",
 //                                symbol="Cs",z=55.,a);  

 //Antimony

 a=121.76*g/mole; 
 G4Element* elSb = new G4Element(name="Antimony",
                                symbol="Sb",z=51.,a);  


//Define Materials
  G4cout << "\nNow Define Materials ..\n" <<G4endl;
  //

//Air at 20 degree C and 1 atm for the ambiet air.
  // Also Air as  a radiator material for inside the tubes.
//--
  density = 1.205e-03*g/cm3;
  pressure=1.*atmosphere;
  temperature=293.*kelvin;
  G4Material* Air = new G4Material(name="Air ", density, numel=2,
                                   kStateGas,temperature,pressure);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  G4double* AirAbsorpLength=new G4double[NumPhotWaveLengthBins];
  G4double* AirRindex=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    AirAbsorpLength[ibin]=1.E32*mm;
    AirRindex[ibin]=1.000273;
  }
  G4MaterialPropertiesTable* AirMPT = 
                            new G4MaterialPropertiesTable();

    AirMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  			AirAbsorpLength,NumPhotWaveLengthBins);

  Air->SetMaterialPropertiesTable(AirMPT);
  RichTbAmbientAir = Air;


  density = 1.205e-03*g/cm3;
  pressure=1.*atmosphere;
  temperature=293.*kelvin;
  G4Material* TAir = new G4Material(name="TAir ", density, numel=2,
                                   kStateGas,temperature,pressure);
  TAir->AddElement(elN, fractionmass=0.7);
  TAir->AddElement(elO, fractionmass=0.3);

  G4double* TAirAbsorpLength=new G4double[NumPhotWaveLengthBins];
  G4double* TAirRindex=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    TAirAbsorpLength[ibin]=1.E32*mm;
    TAirRindex[ibin]=1.000273;
  }
  G4MaterialPropertiesTable* TAirMPT = 
                            new G4MaterialPropertiesTable();

    TAirMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  			TAirAbsorpLength,NumPhotWaveLengthBins);

    TAirMPT->AddProperty("RINDEX", PhotonMomentum, 
                      AirRindex,NumPhotWaveLengthBins);

  TAir->SetMaterialPropertiesTable(TAirMPT);
  RichTbTubeAir = TAir;

  //Nitrogen gas.

  density = 0.8073e-03*g/cm3;
  pressure = RConfig -> getPressureN2();
  temperature = RConfig ->getTemperatureN2();

  G4Material* NitrogenGas = new G4Material(name="NitrogenGas ", 
                                 density, numel=1,
                                 kStateGas,temperature,pressure);
  NitrogenGas->AddElement(elN, natoms=2);

  G4double* NitrogenGasAbsorpLength=new G4double[NumPhotWaveLengthBins];
  G4double* NitrogenGasRindex=new G4double[NumPhotWaveLengthBins];
  G4double* NitrogenGasPhotW=new G4double[NumPhotWaveLengthBins];
  

  std::vector<G4double>N2RefInd= InitN2RefIndex(pressure,temperature);
  std::vector<G4double>N2RefPhotW=InitN2RefPhotW();

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    NitrogenGasAbsorpLength[ibin]=1.E32*mm;

    NitrogenGasRindex[ibin]=N2RefInd[ibin]; 
    NitrogenGasPhotW[ibin]=N2RefPhotW[ibin];

  }
  G4MaterialPropertiesTable* NitrogenGasMPT = 
                            new G4MaterialPropertiesTable();

    NitrogenGasMPT->AddProperty("ABSLENGTH",NitrogenGasPhotW,
  			NitrogenGasAbsorpLength,NumPhotWaveLengthBins);

    NitrogenGasMPT->AddProperty("RINDEX", NitrogenGasPhotW, 
                      NitrogenGasRindex,NumPhotWaveLengthBins);

    NitrogenGas->SetMaterialPropertiesTable(NitrogenGasMPT);
    RichTbNitrogenGas = NitrogenGas;

//Water
 density=1.000*g/cm3;
 G4Material* H2O = new G4Material(name="Water",density,numel=2);
 H2O->AddElement(elH,natoms=2);
 H2O->AddElement(elO,natoms=1);

 G4double* H2OAbsorpLength=new G4double[NumPhotWaveLengthBins];
 G4double* H2ORindex=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    H2OAbsorpLength[ibin]=1.E32*mm;
    H2ORindex[ibin]=1.33;
  }


     G4MaterialPropertiesTable* H2OMPT = 
                           new G4MaterialPropertiesTable();

     H2OMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  		H2OAbsorpLength,NumPhotWaveLengthBins);

    H2OMPT->AddProperty("RINDEX", PhotonMomentum, 
                      H2ORindex,NumPhotWaveLengthBins);

    H2O->SetMaterialPropertiesTable(H2OMPT);


 RichTbH2O=H2O;
//Sio2 
//There is a quartz for the mirror and
 //another quartz which is used in aerogel and
 // yet another quartz used for the quartz window.
 //Mirrorquartz

 density=2.200*g/cm3;
 G4Material* SiO2MirrorQuartz = new G4Material(name="MirrorQuartz",
                                              density,numel=2);
 SiO2MirrorQuartz->AddElement(elSi,natoms=1);
 SiO2MirrorQuartz->AddElement(elO,natoms=2);
 
  G4double* MirrorQuartzRindex=new G4double[NumPhotWaveLengthBins];
  G4double* MirrorQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    MirrorQuartzAbsorpLength[ibin]=0.01*mm;

  }
  G4MaterialPropertiesTable* MirrorQuartzMPT = 
                            new G4MaterialPropertiesTable();


  MirrorQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
 		MirrorQuartzAbsorpLength,NumPhotWaveLengthBins);


  SiO2MirrorQuartz->SetMaterialPropertiesTable(MirrorQuartzMPT);
 RichTbMirrorQuartz=SiO2MirrorQuartz;

 density=2.200*g/cm3;
 G4Material* SiO2AerogelQuartz = new G4Material(name="AerogelQuartz",
                                              density,numel=2);
 SiO2AerogelQuartz->AddElement(elSi,natoms=1);
 SiO2AerogelQuartz->AddElement(elO,natoms=2);

  // QuartzWindow Quartz
 density=2.200*g/cm3;
 G4Material* WindowQuartz = new G4Material(name="WindowQuartz",
                                              density,numel=2);
 WindowQuartz->AddElement(elSi,natoms=1);
 WindowQuartz->AddElement(elO,natoms=2);
 G4double* WindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
 G4double* WindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    WindowQuartzAbsorpLength[ibin]=1.E32*mm;
    WindowQuartzRindex[ibin]=1.4;
  }
  G4MaterialPropertiesTable* WindowQuartzMPT = 
                            new G4MaterialPropertiesTable();

  WindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		WindowQuartzAbsorpLength,NumPhotWaveLengthBins);

  WindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum, 
                        WindowQuartzRindex,NumPhotWaveLengthBins);
  WindowQuartz->SetMaterialPropertiesTable(WindowQuartzMPT);

  RichTbQuartzWindowMaterial=WindowQuartz;
  //for now this is kept to be same as the hpdquartz window.
 density=2.200*g/cm3;
 G4Material* HpdWindowQuartz = new G4Material(name="HpdWindowQuartz",
                                              density,numel=2);
 HpdWindowQuartz->AddElement(elSi,natoms=1);
 HpdWindowQuartz->AddElement(elO,natoms=2);
 G4double* HpdWindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
 G4double* HpdWindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
   HpdWindowQuartzAbsorpLength[ibin]=1.E32*mm;
    HpdWindowQuartzRindex[ibin]=1.40;
  }
  G4MaterialPropertiesTable* HpdWindowQuartzMPT = 
                            new G4MaterialPropertiesTable();

  HpdWindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		HpdWindowQuartzAbsorpLength,NumPhotWaveLengthBins);

  HpdWindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum, 
                        HpdWindowQuartzRindex,NumPhotWaveLengthBins);
  HpdWindowQuartz->SetMaterialPropertiesTable(HpdWindowQuartzMPT);

  HpdQuartzWindowMaterial=HpdWindowQuartz;

  // Borosilcate window of the Pad Hpd
  // for now kept same as the other Hpd.
 density=2.200*g/cm3;
 G4Material* PadHpdWindowQuartz = new G4Material(name="PadHpdWindowQuartz",
                                              density,numel=2);
 PadHpdWindowQuartz->AddElement(elSi,natoms=1);
 PadHpdWindowQuartz->AddElement(elO,natoms=2);
 G4double* PadHpdWindowQuartzRindex=new G4double[NumPhotWaveLengthBins];
 G4double* PadHpdWindowQuartzAbsorpLength=new G4double[NumPhotWaveLengthBins];
 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
   PadHpdWindowQuartzAbsorpLength[ibin]=1.E32*mm;
    PadHpdWindowQuartzRindex[ibin]=1.40;
  }
  G4MaterialPropertiesTable* PadHpdWindowQuartzMPT = 
                            new G4MaterialPropertiesTable();

  PadHpdWindowQuartzMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		PadHpdWindowQuartzAbsorpLength,NumPhotWaveLengthBins);

  PadHpdWindowQuartzMPT->AddProperty("RINDEX", PhotonMomentum, 
                        PadHpdWindowQuartzRindex,NumPhotWaveLengthBins);
  PadHpdWindowQuartz->SetMaterialPropertiesTable(PadHpdWindowQuartzMPT);

  PadHpdQuartzWindowMaterial=PadHpdWindowQuartz;
  //
  //

  G4int filterNumberThisRun=RConfig->GetFilterTNumber();
  // now for the filter material glass d263
 density=2.200*g/cm3;
 G4Material* GlassD263 = new G4Material(name= FilterTypeString[0],
                                              density,numel=2);
 GlassD263->AddElement(elSi,natoms=1);
 GlassD263->AddElement(elO,natoms=2);

 if(filterNumberThisRun >= 0 ) {

   //in the following the +2 is to match the materialproperty bins
   // for the various materials, to avoid the tons of printout from G4.
   // Please see the explanation below for getting the abosorption
   // length of aerogel. The same comments apply here as well.
   // Essentially the measured transmission input here is a combination of
   // the bulk absorption and the fresnel surface loss. One needs to
   // decouple them. Here a partial attempt is made to avoid
   // modifying the G4OpBoundary process. SE. 15-11-2002.
 G4double* GlassD263Rindex=new G4double[NumPhotBinGlassD263Trans+2];
 G4double* GlassD263AbsorpLength=new G4double[NumPhotBinGlassD263Trans+2];
 G4double* GlassD263MomValue = new G4double[NumPhotBinGlassD263Trans+2];
 G4double* currBulkTransFilter = new G4double[NumPhotBinGlassD263Trans+2];
 FilterTrData* CurFil = RConfig->GetFilterTrData();
  std::vector<G4double>GlassD263TransWL = CurFil-> GetTransWL();
  std::vector<G4double>GlassD263Transmis = CurFil->GetTransTotValue();
  G4double FilterHalfZ= CurFil->GetCurFilterThickness();
  
  for (ibin=0; ibin<NumPhotBinGlassD263Trans+2; ibin++){
   
    
    GlassD263Rindex[ibin]=RefIndexGlassD263;
    if(ibin > 0 && ibin < NumPhotBinGlassD263Trans+1 ){
 //now using the formula trans=std::exp(-thickness/absorplength).
      G4int ibina=ibin-1;

    if(GlassD263TransWL[ibina] > 0.0 ) {
    GlassD263MomValue[ibin]= PhotMomWaveConv*eV/GlassD263TransWL[ibina];
    }
     if(GlassD263Transmis[ibina] >0.0 ) {
       // G4double currentfilterRefIndex= GlassD263Rindex[ibin];
       G4double currentAdjacentMediumRefIndex=NitrogenNominalRefIndex;
       // the following needs to be improved in the future 
       // to have a binary search and
       // interpolation between the adjacent  
       // array elements etc. SE 15-11-2002.
       for(size_t ibinr=0; ibinr<N2RefPhotW.size()-1 ; ibinr++){
	 G4double currMomA=GlassD263MomValue[ibin];
         if(currMomA >= N2RefPhotW[ibinr] && currMomA <= N2RefPhotW[ibinr+1]){
	   currentAdjacentMediumRefIndex=N2RefInd[ibinr];
	    }
       
       }

        if( GlassD263Transmis[ibina] > 0.01 ) {   
          currBulkTransFilter[ibin]=
               GetCurrentBulkTrans(GlassD263Rindex[ibin],
       			    currentAdjacentMediumRefIndex,
                            GlassD263Transmis[ibina]);
        } else {
         currBulkTransFilter[ibin]=GlassD263Transmis[ibina];

       }
       if(currBulkTransFilter[ibin] > 0.0 && 
          currBulkTransFilter[ibin] < 0.9995 ) {
        GlassD263AbsorpLength[ibin]=
        -(2.0*FilterHalfZ)/(std::log(currBulkTransFilter[ibin]));
       }else if (currBulkTransFilter[ibin]== 0.0 ) {
         GlassD263AbsorpLength[ibin]=FilterHalfZ/1.0E32;
       }else {
         GlassD263AbsorpLength[ibin]=DBL_MAX;
       }
     }else {

        GlassD263AbsorpLength[ibin]=FilterHalfZ/1.0E32;
     }
    }

  }
      GlassD263MomValue[0]=PhotonMaxEnergy;
      GlassD263AbsorpLength[0]=GlassD263AbsorpLength[1];
      currBulkTransFilter[0]=currBulkTransFilter[1];
   
      G4int mbin=NumPhotBinGlassD263Trans+1;
      GlassD263MomValue[mbin]=PhotonMinEnergy;
      GlassD263AbsorpLength[mbin]=GlassD263AbsorpLength[mbin-1];
      currBulkTransFilter[mbin]=currBulkTransFilter[mbin-1];

  G4MaterialPropertiesTable* GlassD263MPT = 
                            new G4MaterialPropertiesTable();

  GlassD263MPT->AddProperty("ABSLENGTH",GlassD263MomValue,
		GlassD263AbsorpLength,NumPhotBinGlassD263Trans+2);

  GlassD263MPT->AddProperty("RINDEX",GlassD263MomValue, 
                GlassD263Rindex,NumPhotBinGlassD263Trans+2);

  GlassD263->SetMaterialPropertiesTable(GlassD263MPT);
 }

  GlassD263FilterMaterial=GlassD263;
  RichTbFilterMaterial[0]=GlassD263;
  //for the G4Example only 1 filter type is used.
   G4cout << " Now Define Aerogel .." <<G4endl;


//Aerogel upto five types considered so far.
  // in the G4example the same type is repeated 5 times.
//Now for TypeA

  density=0.200*g/cm3;
 
 G4Material* AerogTypeA = 
       new G4Material(name=AerogelTypeString[0], density, numel=2);
 AerogTypeA->AddMaterial(SiO2AerogelQuartz, fractionmass=97.0*perCent);
 AerogTypeA->AddMaterial(H2O, fractionmass=3.0*perCent);


  G4double* AerogTypeARindex=new G4double[NumPhotWaveLengthBins];  
  G4double* AerogTypeAAbsorpLength=new G4double[NumPhotWaveLengthBins];
  G4double* AerogTypeARScatLength = new G4double[NumPhotWaveLengthBins];
  G4double* currentAgelTrans = new G4double[NumPhotWaveLengthBins];

  std::vector<G4double>AerogelTypeASLength = GetAerogelRScatLength(AerogelTypeA);
  G4int AerogNumber=0;
  G4double AerogelLength=GetCurAerogelLength(AerogNumber);
  G4double MaxTotTransmission=AerogelTypeATotTrans;
  // Unfortunately the transmission measurement values only give the 
  // total transmission which includes the loss within aerogel
  // and the Fresnel loss at the surface. In order to 
  // partially decouple this, the approximate loss at the
  // the surface is calculated using the ref index of the
  // aerogel and its surroundings. Then this is added to the
  // measured transmission to get the transmission in the bulk of
  // aerogel. This is then converted to an absorption length. 
  // In a more accurate implementation the loss at the surface
  // should be calculated using a more precise formula. It is 
  // difficult since we do not know the direction of the photons
  // at this point.
  // One possibility is to modify the G4opBoundaryProcess
  // for this, since we do know the direction of the photons by then.
  // This is not done for this G4example, but only in the LHCb implementation.
  // SE 15-11-2002.  
  // The aerogel is inside a volume made of Nitrogen 

 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
   AerogTypeARindex[ibin]= ConvertAgelRIndex(PhotonMomentum[ibin],0);
   AerogTypeARScatLength[ibin]=AerogelTypeASLength[ibin];
   // G4double photwl = PhotMomWaveConv/ (PhotonMomentum[ibin]/eV);

   G4double currentAgelRefIndex=  AerogTypeARindex[ibin];
   G4double currentNeighbourRefIndex= N2RefInd[ibin];
   currentAgelTrans[ibin]= 
        GetCurrentBulkTrans( currentAgelRefIndex,
        currentNeighbourRefIndex,MaxTotTransmission);
 //now using the formula trans=std::exp(-thickness/absorplength)
 // to get the absorplength.

     if( currentAgelTrans[ibin] > 0.0 && currentAgelTrans[ibin] < 0.9995) {
        AerogTypeAAbsorpLength[ibin]=
        -(AerogelLength)/(std::log( currentAgelTrans[ibin]));
     }else if (currentAgelTrans[ibin] == 0.0) {
    
       AerogTypeAAbsorpLength[ibin]=AerogelLength/1.0E32;
     }else {
       
       AerogTypeAAbsorpLength[ibin]=DBL_MAX;
     }

  }

 G4MaterialPropertiesTable* AerogTypeAMPT = 
                          new G4MaterialPropertiesTable();

   AerogTypeAMPT->AddProperty("ABSLENGTH",PhotonMomentum,
 		AerogTypeAAbsorpLength,NumPhotWaveLengthBins);
  
  
  AerogTypeAMPT->AddProperty("RAYLEIGH",PhotonMomentum,
  		     AerogTypeARScatLength,NumPhotWaveLengthBins);

  AerogTypeAMPT->AddProperty("RINDEX", PhotonMomentum, 
                      AerogTypeARindex,NumPhotWaveLengthBins);
   
  AerogTypeA->SetMaterialPropertiesTable(AerogTypeAMPT);


 RichTbAerogelTypeA = AerogTypeA;
 RichTbAerogelMaterial[0] = AerogTypeA;
 // In the G4example the same type is repeated 5 times.
 // in the LHCb implementation 5 types of aerogel materials used.
 //Now for Aerogel TypeB

 RichTbAerogelTypeB = AerogTypeA;
 RichTbAerogelMaterial[1] = AerogTypeA;

 //Now for aerogel TypeC

 RichTbAerogelTypeC = AerogTypeA;
 RichTbAerogelMaterial[2] = AerogTypeA;

 //Now for aerogel TypeD

 RichTbAerogelTypeD = AerogTypeA;
 RichTbAerogelMaterial[3] = AerogTypeA;

 //Now for aerogel Type E


 RichTbAerogelTypeE = AerogTypeA;
 RichTbAerogelMaterial[4] = AerogTypeA;



  //Bialkali Photocathode

 //the following numbers on the property of the BiAlkali Photocathode
  // may not be accurate.
 //Some number is is jut put in for initial program test purposes.
  density=0.100*g/cm3;
 G4Material* BiAlkaliPhCathode = new G4Material(name="BiAlkaliPhCathode", 
                                 density, numel=3);
 BiAlkaliPhCathode->AddElement(elNa, fractionmass=37.5*perCent);
 BiAlkaliPhCathode->AddElement(elK, fractionmass=37.5*perCent);
 BiAlkaliPhCathode->AddElement(elSb, fractionmass=25.0*perCent);

 //for now  properties for the ph cathode material.

 G4double* BiAlkaliPhCathodeRindex=new G4double[NumPhotWaveLengthBins];
 G4double* BiAlkaliPhCathodeAbsorpLength=new G4double[NumPhotWaveLengthBins];
 G4double CathLen=PhotoCathodeThickness;
 G4double CathTrans=PhCathodeNominalTransmission;
 G4double CathAbsorpLen;
 if(CathTrans > 0.0 && CathTrans < 0.9995 ) {
      CathAbsorpLen =  -(CathLen)/(std::log(CathTrans));
 }else if (CathTrans > 0.0) {
     CathAbsorpLen  = CathLen/1.0E32;
 }else {
     CathAbsorpLen  = DBL_MAX; 
 }
 
 for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
   BiAlkaliPhCathodeAbsorpLength[ibin]=CathAbsorpLen;
    BiAlkaliPhCathodeRindex[ibin]=1.40;
  }
  G4MaterialPropertiesTable* BiAlkaliPhCathodeMPT = 
                            new G4MaterialPropertiesTable();
  BiAlkaliPhCathodeMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		BiAlkaliPhCathodeAbsorpLength,NumPhotWaveLengthBins);

  BiAlkaliPhCathodeMPT->AddProperty("RINDEX", PhotonMomentum, 
                        BiAlkaliPhCathodeRindex,NumPhotWaveLengthBins);
  BiAlkaliPhCathode->SetMaterialPropertiesTable(BiAlkaliPhCathodeMPT);
  PadHpdPhCathodeMaterial=BiAlkaliPhCathode;

//CF4
//no data available at room temp and pressure;
 density=0.003884*g/cm3;
 temperature=273.*kelvin;
 pressure=1.0*atmosphere;
 a=88.01*g/mole;

 G4Material* CF4 =new G4Material(name="CF4",density,numel=2,
                             kStateGas,temperature,pressure);
 CF4->AddElement(elC,natoms=1);
 CF4->AddElement(elF,natoms=4);
 // Sellmeir coef to be added.
 RichTbCF4=CF4;

  G4cout << "\nNowDefineVacuum ..\n" <<G4endl;

//Vacuum
//
 density=universe_mean_density;   
 a=1.01*g/mole;
 pressure=1.e-19*pascal;
 temperature=0.1*kelvin;

 G4Material* vacuum = new G4Material(name="Galactic",density,numel=1,
                              kStateGas,temperature,pressure);
 vacuum->AddElement(elH,natoms=1);

  G4double* VacAbsorpLength=new G4double[NumPhotWaveLengthBins];
  G4double* VacRindex=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    VacAbsorpLength[ibin]=1.E32*mm;
    // the following ref index is just artifical, just to
    // avoid the refraction between nitrogen gas and hpd master.
    VacRindex[ibin]=1.000273;
  }
  G4MaterialPropertiesTable* VacMPT = 
                            new G4MaterialPropertiesTable();

  VacMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		      VacAbsorpLength,NumPhotWaveLengthBins);
  VacMPT->AddProperty("RINDEX", PhotonMomentum, 
                      VacRindex,NumPhotWaveLengthBins);
  vacuum->SetMaterialPropertiesTable(VacMPT);
    
  RichTbVacuum=vacuum;

//beamgas
//
  density=1.e-5*g/cm3;
  pressure=2.e-2*bar;
  temperature=STP_Temperature;   
  G4Material* beamgas = new G4Material(name="Beamgas",density,numel=1,
          		         kStateGas,temperature,pressure);
  beamgas->AddMaterial(Air,fractionmass=1.); // beware that air is at 20 deg;
  
//
//Aluminium
  density=2.7*g/cm3;
  G4Material* Aluminium =new G4Material(name="Aluminium",density,numel=1);
  Aluminium->AddElement(elAL,natoms=1);

  G4double* AluminiumAbsorpLength=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    AluminiumAbsorpLength[ibin]=0.0*mm;
  }

  G4MaterialPropertiesTable* AluminiumMPT = 
                            new G4MaterialPropertiesTable();

  AluminiumMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		AluminiumAbsorpLength,NumPhotWaveLengthBins);

  Aluminium->SetMaterialPropertiesTable(AluminiumMPT);
  RichTbAluminium=Aluminium;
//PlasticAg , this is used as a wrap of aerogel and as upstream holder
  // for aerogel frame. For now use same properties as that of Aluminium.
  // this is just an opaque material.

  density=2.7*g/cm3;
  G4Material* PlasticAg =new G4Material(name="PlasticAg",density,numel=1);
  PlasticAg->AddElement(elAL,natoms=1);

  G4double* PlasticAgAbsorpLength=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    PlasticAgAbsorpLength[ibin]=0.0*mm;
  }

  G4MaterialPropertiesTable* PlasticAgMPT = 
                            new G4MaterialPropertiesTable();

  PlasticAgMPT->AddProperty("ABSLENGTH",PhotonMomentum,
	        PlasticAgAbsorpLength,NumPhotWaveLengthBins);

  PlasticAg->SetMaterialPropertiesTable(PlasticAgMPT);
  RichTbPlasticAg=PlasticAg;
// Kovar
  density=2.7*g/cm3;
  G4Material* Kovar =new G4Material(name="Kovar",density,numel=1);
  Kovar->AddElement(elAL,natoms=1);

  G4double* KovarAbsorpLength=new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    KovarAbsorpLength[ibin]=0.0*mm;
  }

  G4MaterialPropertiesTable* KovarMPT = 
                            new G4MaterialPropertiesTable();

  KovarMPT->AddProperty("ABSLENGTH",PhotonMomentum,
		KovarAbsorpLength,NumPhotWaveLengthBins);

  Kovar->SetMaterialPropertiesTable(KovarMPT);
  HpdTubeMaterial=Kovar;

// Silicon

  density=2.33*g/cm3;
  G4Material* Silicon =new G4Material(name="Silicon",density,numel=1);
  Silicon->AddElement(elSi,natoms=1);

  G4double* SiliconAbsorpLength=new G4double[NumPhotWaveLengthBins];

   for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    SiliconAbsorpLength[ibin]=0.0*mm;
    }

  G4MaterialPropertiesTable* SiliconMPT = 
                            new G4MaterialPropertiesTable();

   SiliconMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  		SiliconAbsorpLength,NumPhotWaveLengthBins);

  Silicon->SetMaterialPropertiesTable(SiliconMPT);
  HpdSiDetMaterial=Silicon;

// Silicon coating made of Si02.

  density=2.33*g/cm3;
  G4Material* SiliconCoating =new G4Material(name="SilCoat",density,numel=2);
  SiliconCoating->AddElement(elSi,natoms=1);
  SiliconCoating->AddElement(elO,natoms=2);

  G4double* SiliconCoatingAbsorpLength=new G4double[NumPhotWaveLengthBins];
  // G4double* SiliconCoatingRindex=new G4double[NumPhotWaveLengthBins];

   for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    SiliconCoatingAbsorpLength[ibin]=0.0001*mm;

    }

  G4MaterialPropertiesTable* SiliconCoatingMPT = 
                            new G4MaterialPropertiesTable();

   SiliconCoatingMPT->AddProperty("ABSLENGTH",PhotonMomentum,
  		SiliconCoatingAbsorpLength,NumPhotWaveLengthBins);

  SiliconCoating->SetMaterialPropertiesTable(SiliconCoatingMPT);
  HpdSiCoatingMaterial=SiliconCoating;

  //
  // Now for the material properties of Surfaces
  //
  //
  //
  //Front (reflecting surface of RichTb Mirror)

  // First define wavelength in nm.
  //For now assume that all segments have the same reflectivity.
  // Hence the reflectivity is defined outside the loop of the
  // the number of segments.
  //Only the front surface is created.
  // The abosorption length is set to a small value just to
  // avoid photons exiting from the back of the mirror.
  // the efficiency is for the absorption process.

  
  G4double* PhotonMomentumRefl
    =new G4double[NumPhotonRichMirrorReflWaveLengthBins];
  G4double* PhotWaveRefl = 
       new  G4double[NumPhotonRichMirrorReflWaveLengthBins];
  G4double* PhotReflEff =new  G4double[NumPhotonRichMirrorReflWaveLengthBins]; 
  G4double* MirrorQuRefIndex 
         =new  G4double[NumPhotonRichMirrorReflWaveLengthBins]; 

  for (ibin=0; ibin<NumPhotonRichMirrorReflWaveLengthBins; ibin++){
    PhotonMomentumRefl[ibin]=PhotMomWaveConv*eV/ PhotonWavelengthRefl[ibin];
     PhotWaveRefl[ibin]=  RichTbMirrorReflectivity[ibin];
     PhotReflEff[ibin]= RichTbMirrorEfficiency[ibin];
    //the following lines to avoid reflection at the mirror.

    MirrorQuRefIndex[ibin] = 1.40;
  }

   G4OpticalSurface * OpRichTbMirrorSurface =
     new G4OpticalSurface("RichTbMirrorSurface");

   OpRichTbMirrorSurface->SetType(dielectric_metal);
   OpRichTbMirrorSurface->SetFinish(polished);
   OpRichTbMirrorSurface->SetModel(glisur);
   G4MaterialPropertiesTable* OpRichTbMirrorSurfaceMPT = 
                            new G4MaterialPropertiesTable();

       OpRichTbMirrorSurfaceMPT->AddProperty("REFLECTIVITY",
                          PhotonMomentumRefl,
    			   PhotWaveRefl,
                          NumPhotonRichMirrorReflWaveLengthBins);
     OpRichTbMirrorSurfaceMPT->AddProperty("EFFICIENCY",
                          PhotonMomentumRefl,
    			   PhotReflEff,
                          NumPhotonRichMirrorReflWaveLengthBins);
    OpRichTbMirrorSurfaceMPT->AddProperty("RINDEX",
                          PhotonMomentumRefl,
  			  MirrorQuRefIndex,
                          NumPhotonRichMirrorReflWaveLengthBins);

  OpRichTbMirrorSurface->SetMaterialPropertiesTable(OpRichTbMirrorSurfaceMPT);
  RichTbOpticalMirrorSurface=OpRichTbMirrorSurface;


  //  OpRichTbMirrorSurface->DumpInfo();

  // Now for the Surface of the Vessel Enclosure.


   G4OpticalSurface * OpRichTbEnclosureSurface =
     new G4OpticalSurface("RichTbEnclosureSurface");
   OpRichTbEnclosureSurface->SetType(dielectric_metal);
   OpRichTbEnclosureSurface->SetFinish(polished);
   OpRichTbEnclosureSurface->SetModel(glisur);

   G4double NumPhotonRichEnclosureSurfaceWaveLengthBins=10;
   G4double  RichTbEnclosureSurfaceReflectivity[]=
   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double RichTbEnclosureSurfaceEfficiency[]=
   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   G4double RichEnclosureSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
 
   G4MaterialPropertiesTable* OpRichTbEnclosureSurfaceMPT = 
                            new G4MaterialPropertiesTable();

   OpRichTbEnclosureSurfaceMPT->AddProperty("REFLECTIVITY",
                            RichEnclosureSurfacePhotMom,
			   RichTbEnclosureSurfaceReflectivity,
                           static_cast<int>(NumPhotonRichEnclosureSurfaceWaveLengthBins));
   OpRichTbEnclosureSurfaceMPT->AddProperty("EFFICIENCY",
                           RichEnclosureSurfacePhotMom,
			   RichTbEnclosureSurfaceEfficiency,
                           static_cast<int>(NumPhotonRichEnclosureSurfaceWaveLengthBins));

  OpRichTbEnclosureSurface->
        SetMaterialPropertiesTable(OpRichTbEnclosureSurfaceMPT);

  RichTbOpticalEnclosureSurface=OpRichTbEnclosureSurface;

  //Now for the surface between the TAir and Quartz Window of the HPD

   G4OpticalSurface * OpHpdQuartzWTSurface =
     new G4OpticalSurface("HpdQuartzWTSurface");
   OpHpdQuartzWTSurface->SetType(dielectric_dielectric);
   OpHpdQuartzWTSurface->SetFinish(polished);
   OpHpdQuartzWTSurface->SetModel(glisur);
   //OpHpdQuartzWTSurface->SetModel(unified);

   G4double NumPhotonHpdQuartzWTSurfaceWaveLengthBins=10;
   G4double  HpdQuartzWTSurfaceReflectivity[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double HpdQuartzWTSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double HpdQuartzWTSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
 
   G4MaterialPropertiesTable* OpHpdQuartzWTSurfaceMPT = 
                            new G4MaterialPropertiesTable();


   OpHpdQuartzWTSurfaceMPT->AddProperty("REFLECTIVITY",
                            HpdQuartzWTSurfacePhotMom,
			   HpdQuartzWTSurfaceReflectivity,
                           static_cast<int>(NumPhotonHpdQuartzWTSurfaceWaveLengthBins));
   OpHpdQuartzWTSurfaceMPT->AddProperty("EFFICIENCY",
                           HpdQuartzWTSurfacePhotMom,
			   HpdQuartzWTSurfaceEfficiency,
                           static_cast<int>(NumPhotonHpdQuartzWTSurfaceWaveLengthBins));

  OpHpdQuartzWTSurface->
        SetMaterialPropertiesTable(OpHpdQuartzWTSurfaceMPT);

   HpdTQuartzWSurface=OpHpdQuartzWTSurface;



  //Now for the surface between the Quartz Window and Ph cathode of the HPD

   G4OpticalSurface * OpHpdQuartzWPSurface =
     new G4OpticalSurface("HpdQuartzWPSurface");
   OpHpdQuartzWPSurface->SetType(dielectric_dielectric);
   OpHpdQuartzWPSurface->SetFinish(polished);
   OpHpdQuartzWPSurface->SetModel(glisur);

   G4double NumPhotonHpdQuartzWPSurfaceWaveLengthBins=10;
   G4double  HpdQuartzWPSurfaceReflectivity[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


   G4double HpdQuartzWPSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double HpdQuartzWPSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
 
   G4MaterialPropertiesTable* OpHpdQuartzWPSurfaceMPT = 
                            new G4MaterialPropertiesTable();


   OpHpdQuartzWPSurfaceMPT->AddProperty("REFLECTIVITY",
                            HpdQuartzWPSurfacePhotMom,
			   HpdQuartzWPSurfaceReflectivity,
                           static_cast<int>(NumPhotonHpdQuartzWPSurfaceWaveLengthBins));
   OpHpdQuartzWPSurfaceMPT->AddProperty("EFFICIENCY",
                           HpdQuartzWPSurfacePhotMom,
			   HpdQuartzWPSurfaceEfficiency,
                           static_cast<int>(NumPhotonHpdQuartzWPSurfaceWaveLengthBins));

  OpHpdQuartzWPSurface->
        SetMaterialPropertiesTable(OpHpdQuartzWPSurfaceMPT);

   HpdQuartzWPhCathodeSurface=OpHpdQuartzWPSurface;



  //Now for the skin surface of the PhCathode so that photons do
  // not come out of the Photocathode.
   // Changed to dielectric-dielectric so that photons DO come out
   // of the photocathode. SE 26-9-01.

   G4OpticalSurface * OpPhCathodeSurface =
     new G4OpticalSurface("PhCathodeSurface");

   OpPhCathodeSurface->SetType(dielectric_dielectric);
   OpPhCathodeSurface->SetFinish(polished);
   OpPhCathodeSurface->SetModel(glisur);


   G4double NumPhotonPhCathodeSurfaceWaveLengthBins=10;
   G4double  PhCathodeSurfaceReflectivity[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


   G4double PhCathodeSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double PhCathodeSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
 
   G4MaterialPropertiesTable* OpPhCathodeSurfaceMPT = 
                            new G4MaterialPropertiesTable();


   OpPhCathodeSurfaceMPT->AddProperty("REFLECTIVITY",
                           PhCathodeSurfacePhotMom,
			   PhCathodeSurfaceReflectivity,
                           static_cast<int>(NumPhotonPhCathodeSurfaceWaveLengthBins));
   OpPhCathodeSurfaceMPT->AddProperty("EFFICIENCY",
                           PhCathodeSurfacePhotMom,
			   PhCathodeSurfaceEfficiency,
                           static_cast<int>(NumPhotonPhCathodeSurfaceWaveLengthBins));

  OpPhCathodeSurface->
        SetMaterialPropertiesTable(OpPhCathodeSurfaceMPT);

  PhCathodeSkinSurface=OpPhCathodeSurface;
  PhCathodeBorderSurface=OpPhCathodeSurface;




  //Now for the surface between Interior of HPD and Silicon Coating.

   G4OpticalSurface * OpHpdSiCoatSurface =
     new G4OpticalSurface("HpdSiCoatSurface");
   OpHpdSiCoatSurface->SetType(dielectric_metal);
   OpHpdSiCoatSurface->SetFinish(polished);
   OpHpdSiCoatSurface->SetModel(glisur);


   G4double NumPhotonHpdSiCoatSurfaceWaveLengthBins=10;

  G4double  HpdSiCoatSurfaceReflectivity[]=
    {0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9};

   G4double HpdSiCoatSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double HpdSiCoatSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
   G4double HpdSiCoatSurfaceRefInd[]=
    {1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4};
 
   
   G4MaterialPropertiesTable* OpHpdSiCoatSurfaceMPT = 
                            new G4MaterialPropertiesTable();


   OpHpdSiCoatSurfaceMPT->AddProperty("REFLECTIVITY",
                            HpdSiCoatSurfacePhotMom,
			   HpdSiCoatSurfaceReflectivity,
                           static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));
   OpHpdSiCoatSurfaceMPT->AddProperty("EFFICIENCY",
                           HpdSiCoatSurfacePhotMom,
			   HpdSiCoatSurfaceEfficiency,
                           static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));
   OpHpdSiCoatSurfaceMPT->AddProperty("RINDEX",
                           HpdSiCoatSurfacePhotMom,
			   HpdSiCoatSurfaceRefInd,
                           static_cast<int>(NumPhotonHpdSiCoatSurfaceWaveLengthBins));

  OpHpdSiCoatSurface->
        SetMaterialPropertiesTable(OpHpdSiCoatSurfaceMPT);

   HpdSiCoatSurface=OpHpdSiCoatSurface;



  // Now for the Surface of the MetalTube of HPD.


   G4OpticalSurface * OpRichTbHpdMetalSurface =
     new G4OpticalSurface("RichTbHpdMetalSurface");
   OpRichTbHpdMetalSurface->SetType(dielectric_metal);
   OpRichTbHpdMetalSurface->SetFinish(polished);
   OpRichTbHpdMetalSurface->SetModel(glisur);

   G4double NumPhotonHpdMetalSurfaceWaveLengthBins=10;
   G4double  RichHpdMetalSurfaceReflectivity[]=
   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

   G4double RichHpdMetalSurfaceEfficiency[]=
   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   G4double RichHpdMetalSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
 
   G4MaterialPropertiesTable* OpRichTbHpdMetalSurfaceMPT = 
                            new G4MaterialPropertiesTable();

   OpRichTbHpdMetalSurfaceMPT->AddProperty("REFLECTIVITY",
                            RichHpdMetalSurfacePhotMom,
			   RichHpdMetalSurfaceReflectivity,
                           static_cast<int>(NumPhotonHpdMetalSurfaceWaveLengthBins));
   OpRichTbHpdMetalSurfaceMPT->AddProperty("EFFICIENCY",
                           RichHpdMetalSurfacePhotMom,
			   RichHpdMetalSurfaceEfficiency,
                           static_cast<int>(NumPhotonHpdMetalSurfaceWaveLengthBins));

  OpRichTbHpdMetalSurface->
        SetMaterialPropertiesTable(OpRichTbHpdMetalSurfaceMPT);

  RichTbOpticalHpdMetalSurface=OpRichTbHpdMetalSurface;


  //Now for the surface of the Filter

   G4OpticalSurface * OpRichTbFilterSurface =
     new G4OpticalSurface("RichTbFilterSurface");
     OpRichTbFilterSurface->SetType(dielectric_dielectric);
     OpRichTbFilterSurface->SetFinish(polished);
     OpRichTbFilterSurface->SetModel(glisur);



 if(filterNumberThisRun >= 0 ) {

   G4int FilterNumbins=NumPhotonRichTbFilterSurfaceWaveLengthBins;


   G4double*  FilterReflectivity = new G4double(FilterNumbins);
   G4double* FilterEff =new G4double(FilterNumbins);
   G4double* FilterPhotMom =new G4double(FilterNumbins);


     for(G4int ibinf =0 ; ibinf < FilterNumbins; ibinf++ ){
        FilterReflectivity[ibinf]= RichTbFilterSurfaceReflectivity[ibinf];
        FilterEff[ibinf]= RichTbFilterSurfaceEfficiency[ibinf];
        FilterPhotMom[ibinf]= RichTbFilterSurfacePhotMom[ibinf];

//       G4MaterialPropertiesTable* OpRichTbFilterSurfaceMPT = 
//                         new G4MaterialPropertiesTable();

     }
   RichTbOpticalFilterSurface=OpRichTbFilterSurface;

 }

  delete [] PhotonMomentum;
  delete [] AirAbsorpLength;
  delete [] AirRindex;
  delete [] MirrorQuartzRindex;
  delete [] MirrorQuartzAbsorpLength;
  delete [] WindowQuartzRindex;
  delete [] WindowQuartzAbsorpLength;
  delete [] AluminiumAbsorpLength;
  delete [] KovarAbsorpLength;
  delete [] PhotonMomentumRefl;


}
G4double RichTbMaterial::ConvertAgelRIndex(G4double phmom, G4int AgelTnum ) {
  AerogelRefData* AgData= rConfig -> GetAerogelRefdata();
  //Now to convert and interpolate to get the same binning
  // as the other property vectors.
  G4double Refind=0.;
  G4int Numphbin=AgData-> GetNumberOfRefIndBins();
  G4double phm1,phm2;
  if(phmom < AgData->GetAerogelRefphotE(0) ){
    Refind=AgData->GetCurAerogelRefIndValue(0,AgelTnum ); }
  if(phmom >= AgData->GetAerogelRefphotE(Numphbin-1 ) ) {
  Refind=AgData->GetCurAerogelRefIndValue(Numphbin-1,AgelTnum ); }

  for( G4int iba=0; iba<Numphbin-1 ; iba ++ ) {
   
    phm1=AgData->GetAerogelRefphotE(iba);
    phm2=AgData->GetAerogelRefphotE(iba+1);
 
    if(phmom >= phm1 && phmom < phm2 ) {
      
      G4double ref1=AgData->GetCurAerogelRefIndValue(iba,AgelTnum );
      G4double ref2=AgData->GetCurAerogelRefIndValue(iba+1,AgelTnum );

      G4double grad = (ref2-ref1)/(phm2-phm1);
      G4double aint = ref1- grad*phm1;
      Refind = grad*phmom + aint ;
      break;
    }   
  }

  return Refind;
}
RichTbMaterial::RichTbMaterial() { ; }
RichTbMaterial::~RichTbMaterial(){ ; }



