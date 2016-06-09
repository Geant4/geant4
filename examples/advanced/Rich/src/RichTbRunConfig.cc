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
// RichTbRunConfig.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "RichTbRunConfig.hh"
#include "RichTbMaterialParameters.hh"
#include "globals.hh"
RichTbRunConfig::RichTbRunConfig()
 :CurAerogelTNumber(std::vector<G4int>(MaxNumberOfAerogelTiles)),
  CurAerogelType(std::vector<AerogelType>(MaxNumberOfAerogelTiles)){ 

  const char* RunConfigFile = "RunConfig.dat" ;
 G4cout<<" Run Configuration Input is from "<< RunConfigFile<<G4endl; 
 G4double rPressureN2,rTemperatureN2;
 G4double rMirrorAddTiltXDeg;
 G4double rMirrorAddTiltYDeg;
 G4double unitofPressureN2=1.0*bar;
 G4double unitofTemperatureN2=1.0*kelvin;

 // To reduce the number of control variables to input 
 // from the RunConfig.dat
 // file, define a defualt set for most them here for the G4Example.
  // Pressure in bar for the nitrogen.
 rPressureN2=1.0;
 //temperature of nitrogen in kelvin.
 rTemperatureN2=292.0;
 // refractive index parametrization for nitrogen.
 // to be removed.
 RefIndexPara=1;
 // Mirror tilt in X and Y
 rMirrorAddTiltXDeg=0.0;
 rMirrorAddTiltYDeg=0.0;

 // Assign a  default for the graphics variables.
 RichTbHall_visib=0;
 RichTbEnclosure_visib=2;
 RichTbRadFrame_visib=2;
 RichTbRadUpW_visib=1;
 RichTbRadDnW_visib=1;
 RichTbAerogel_visib=1;
 RichTbAerogelWrap_visib=1;
 RichTbFilter_visib=1;
 RichTbMirror_visib=1;
 RichTbHpdMaster_visib=1;
 RichTbHpdEnvelopeTube_visib=0;
 RichTbHpdQuartzW_visib=1;
 RichTbHpdPhCathode_visib=1;
 RichTbHpdSiDet_visib=1;
 RichTbHpdSectCoat_visib=0;
 RichTbHpdSiPx_visib=0;
 RichTbDrawTrajectory_visib=1;
 // Photoelectron energy in keV
 HpdPhElectronEnergy=16.0;
 FilterTNumber=1;
 NumberOfAerogelTiles=1;
 CurAerogelTNumber[0]=0;
 //
 RichTbNumPartEvent=1;
 RichTbParticleTypeCode=0;
 RichTbParticleStartPosCode=0;
 RichTbParticleDirectionCode=0;
 RichTbParticleEnergyCode=0;
 RichTbParticleEnergy=9.0;
 RichTbPhotLowE=1.305;
 RichTbPhotHighE=7.2;
 WriteMiscOutputFlag=0;
 OutputFileName ="../outputData/MyOutputDataFile";
 OutputHistoDirName="../HistoOutputVop/";
 AerogelRefDataInputFileName="../inputData/aerogelRefIndexInput.dat";
 FilterTransDataInputFileName="../inputData/FilterTransInput.dat";
 MiscOutFileName="../outputData/MyMiscOutputFile";
//
  if(RunConfigFile == 0 ) { 
    G4cout<<"Unable to read the run configuration file "
          <<" using default values for the control variables "
          <<" and input and output file names "<<G4endl; 
  }else {
  G4cout<<"Now getting the run configuration file. "
        << "the Runconfig file name is "<< RunConfigFile <<G4endl;

  std::ifstream finp(RunConfigFile);

 finp>>RichTbHall_visib;
 finp>>RichTbEnclosure_visib;
 finp>>RichTbRadFrame_visib;
 finp>>RichTbRadUpW_visib;
 finp>>RichTbRadDnW_visib;
 finp>>RichTbAerogel_visib;
 finp>>RichTbAerogelWrap_visib;
 finp>>RichTbFilter_visib;
 finp>>RichTbMirror_visib;
 finp>>RichTbHpdMaster_visib;
 finp>>RichTbHpdEnvelopeTube_visib;
 finp>>RichTbHpdQuartzW_visib;
 finp>>RichTbHpdPhCathode_visib;
 finp>>RichTbHpdSiDet_visib;
 finp>>RichTbHpdSectCoat_visib;
 finp>>RichTbHpdSiPx_visib;
 finp>>RichTbDrawTrajectory_visib; 
 finp>>RichTbParticleEnergy;
 finp>>FilterTNumber;
 finp>>WriteOutputFile;
 finp>>OutputFileName;
 finp>>OutputHistoDirName;
 finp>>AerogelRefDataInputFileName;
 finp>>FilterTransDataInputFileName;
  }
 CurrentPressureN2=rPressureN2*unitofPressureN2;
 CurrentTemperatureN2=rTemperatureN2*unitofTemperatureN2;
 MirrorAddTiltX = (rMirrorAddTiltXDeg*pi/180)*rad;
 MirrorAddTiltY = (rMirrorAddTiltYDeg*pi/180)*rad;
 G4cout<<"Current Run Configuration "<<G4endl;
 G4cout<<"Nitrogen Gas  Pressure=  "<<
    CurrentPressureN2/unitofPressureN2<<"  bar"<<G4endl;
 G4cout<<"Nitrogen Gas  Temperature =  "
       <<CurrentTemperatureN2/unitofTemperatureN2<<
        "   kelvin"<<G4endl;
 G4cout<<"Ref Index parametrization used "<<RefIndexPara<<G4endl;
 G4cout<<"Addtional Mirror tilts wrt X and Y in rad "<<MirrorAddTiltX 
       <<"  "<<MirrorAddTiltY<<G4endl;
  //For Following variables 0 means make the volume invisible;
  //                        1 means make it visible as a solid.
  //                        2 means make it visible as a wireframe.
 G4cout<<" RichTbHall_visib = "<<RichTbHall_visib<<G4endl;
 G4cout<<"RichTbEnclosure_visib = "<<RichTbEnclosure_visib<<G4endl;
 G4cout<<"RichTbRad Frame_visib, FrameUpW_visib, FrameDnW  = "
        <<RichTbRadFrame_visib
       <<"   "<<RichTbRadUpW_visib<<"  "<<RichTbRadDnW_visib
       <<G4endl;
 G4cout<<"RichTbAerogel_visib, AerogelWrap_visib  = "
       <<RichTbAerogel_visib<<RichTbAerogelWrap_visib<<G4endl;
 G4cout<<"RichTbFilter_visib = "<<RichTbFilter_visib<<G4endl;
 G4cout<<"RichTbMirror_visib = "<<RichTbMirror_visib<<G4endl;
 G4cout<<"RichTbHpdMaster_visib = "<<RichTbHpdMaster_visib<<G4endl;
 G4cout<<"RichTbHpdEnvelopeTube_visib = "<<RichTbHpdEnvelopeTube_visib<<G4endl;
 G4cout<<"RichTbHpdQuartzW_visib = "<<RichTbHpdQuartzW_visib<<G4endl;
 G4cout<<"RichTbHpdPhCathode_visib = "<<RichTbHpdPhCathode_visib<<G4endl;
 G4cout<<"RichTbHpdSiDet_visib = "<<RichTbHpdSiDet_visib<<G4endl;
 G4cout<<"RichTbHpdSiSectCoat_visib = "<<RichTbHpdSectCoat_visib<<G4endl;
 G4cout<<"RichTbHpdSiPx_visib = "<<RichTbHpdSiPx_visib<<G4endl;
 G4cout<<"RichTbDrawTrajectory_visib =  "<<RichTbDrawTrajectory_visib<<G4endl;
 G4cout<<"RichTbPhElectronEnergy = "<< HpdPhElectronEnergy<<G4endl;

 // In the G4example only upto 1 type of filter is used.
 if (FilterTNumber > 0 ) {
   G4cout<<"RichTbRunConfig: In the G4example "
         <<" only 1 type of filter used " <<G4endl;
   FilterTNumber =0;
 }


 FilterTransData = 
   new FilterTrData(FilterTNumber,FilterTransDataInputFileName);
 CurFilterType=  FilterTransData->GetFilterTypeIndex();

 G4cout<<"Filter Number   =  "<<FilterTNumber<<G4endl;
 if(FilterTNumber >= 0 ) {
   G4cout<<"  Filter Type = "<<CurFilterType<<G4endl;
 } else {
   G4cout<<"No Filter for this Run "<<G4endl;
 }

 // for the G4example 1 tile is used. In the LHCb implementation
 // upto 5 tiles can be used. 
 if(NumberOfAerogelTiles > 1 ) {
   G4cout<<" RichTbRunConfig: For the G4Example " 
         <<" the number of aerogel tiles is limited to 1  "
         <<" Hence setting the number of Tiles to 1 "<<G4endl;
   NumberOfAerogelTiles =1;
 }
 for(G4int iia=0; iia<NumberOfAerogelTiles; iia++) {
   if(CurAerogelTNumber[iia] == 0 ) CurAerogelType[iia]=AerogelTypeA; 
   if(CurAerogelTNumber[iia] == 1 ) CurAerogelType[iia]=AerogelTypeB; 
   if(CurAerogelTNumber[iia] == 2 ) CurAerogelType[iia]=AerogelTypeC; 
   if(CurAerogelTNumber[iia] == 3 ) CurAerogelType[iia]=AerogelTypeD; 
   if(CurAerogelTNumber[iia] == 4 ) CurAerogelType[iia]=AerogelTypeE; 
    }
 G4cout<<"Number of Aerogel Tiles "<<NumberOfAerogelTiles<<G4endl;
 for(G4int ia=0; ia<NumberOfAerogelTiles; ia++ ) {
   G4cout<<" AerogelNumber and Type "<<CurAerogelTNumber[ia]<<"   "<<
     CurAerogelType[ia]<<G4endl;
   }

 G4cout<<"Aerogel RefIndex Input file Name is "
       <<AerogelRefDataInputFileName<<G4endl;
 
 AerogelRefractiveInd= new AerogelRefData(AerogelRefDataInputFileName);

 G4cout<<"Write Outfile File Flag set to "<<WriteOutputFile<<G4endl;
 G4cout<<"OutputData file = "<<OutputFileName <<G4endl;
 G4cout<<"OutputHistoDir = "<<OutputHistoDirName <<G4endl;

 G4String PaType;
 if(RichTbParticleTypeCode == 0 )PaType="piminus";
 if(RichTbParticleTypeCode == 1 )PaType="optical photon";
 if(RichTbParticleTypeCode == 2 )PaType="Mixture of piminus and proton";

 G4cout<<"Number of beam particles per event = "<<RichTbNumPartEvent<<G4endl;

 G4cout<<"Particle Type Code "<<RichTbParticleTypeCode
       <<"  "<< PaType <<G4endl;

 G4cout<<"Particle Start Pos and Direction Energy codes "
   << RichTbParticleStartPosCode <<"  "<< RichTbParticleDirectionCode 
       <<"   "<<RichTbParticleEnergyCode<<G4endl;
 G4cout<<"Energy code 0 is for charged particle and 1 is for photon generation"
       << " using particle gun "<<G4endl;
 G4cout<<"Gun Gen Charged Particle energy in GeV is "
       <<RichTbParticleEnergy<<G4endl;
 G4cout<<"Gun Gen Photon energy range in eV is "<<RichTbPhotLowE
       <<"  "<< RichTbPhotHighE <<G4endl;

 }
RichTbRunConfig::~RichTbRunConfig() { ; }






