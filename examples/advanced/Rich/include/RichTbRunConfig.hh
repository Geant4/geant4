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
// RichTbRunConfig.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbRunConfig_h
#define RichTbRunConfig_h 1

#include "globals.hh"
#include <vector>
#include "AerogelTypeSpec.hh"
#include "FilterTypeSpec.hh"
#include "AerogelRefData.hh"
#include "FilterTrData.hh"
#include <cmath>
#include <fstream>
class RichTbRunConfig{
 public:
  RichTbRunConfig();
 virtual ~RichTbRunConfig();
  G4double getPressureN2()
  {return CurrentPressureN2;}
  G4double getTemperatureN2()
   { return CurrentTemperatureN2; }
  G4double getMirrorAddTiltX()
  {return MirrorAddTiltX; }
  G4double getMirrorAddTiltY()
  {return MirrorAddTiltY; }
  G4int getRichTbHall_visib()
  {return  RichTbHall_visib;}
  G4int getRichTbEnclosure_visib()
  {return  RichTbEnclosure_visib;}
  G4int getRichTbRadFrame_visib()
  {return  RichTbRadFrame_visib;}
  G4int getRichTbRadUpW_visib()
  {return  RichTbRadUpW_visib;}
  G4int getRichTbRadDnW_visib()
  {return  RichTbRadDnW_visib;}
  G4int getRichTbAerogel_visib()
  {return  RichTbAerogel_visib;}
  G4int getRichTbAerogelWrap_visib()
  {return  RichTbAerogelWrap_visib;}
  G4int getRichTbFilter_visib()
  {return  RichTbFilter_visib;}
  G4int getRichTbMirror_visib()
  {return RichTbMirror_visib;}
  G4int getRichTbHpdMaster_visib()
  {return  RichTbHpdMaster_visib; }
  G4int  getRichTbHpdEnvelopeTube_visib()
  {return  RichTbHpdEnvelopeTube_visib;}

  G4int getRichTbHpdQuartzW_visib() 
  {return RichTbHpdQuartzW_visib; }
  G4int getRichTbHpdPhCathode_visib() 
  {return RichTbHpdPhCathode_visib; }
  G4int getRichTbHpdSiDet_visib() 
  {return RichTbHpdSiDet_visib; }
  G4int getRichTbHpdSectCoat_visib() 
  {return RichTbHpdSectCoat_visib; }
  G4int getRichTbHpdSiPx_visib() 
  {return RichTbHpdSiPx_visib; }
  G4int getRefIndexPara() 
  {return RefIndexPara; }
  G4int getRichTbDrawTrajectory()  {return RichTbDrawTrajectory_visib; }
  G4double getHpdPhElectronEnergy() {return HpdPhElectronEnergy; }
  G4int  DoWriteOutputFile() {return WriteOutputFile; }
  G4String getOutputFileName() {return OutputFileName; }
  G4String getOutputHistoDirName() {return OutputHistoDirName; }
  G4int getWriteMiscOutputFlag() {return  WriteMiscOutputFlag; }
  G4String getMiscOutFileName() {return MiscOutFileName; }
  G4int GetFilterTNumber() {return FilterTNumber; }
  FilterType GetFilterType() {return CurFilterType; }
  G4int GetNumberOfAerogelTiles() {return NumberOfAerogelTiles ; }
  G4int GetCurAerogelTNumber(G4int AgelNum ) 
    {return CurAerogelTNumber[AgelNum];}
  AerogelType GetCurAerogelType(G4int AgelNumber) 
         {return CurAerogelType[AgelNumber];}
  AerogelRefData* GetAerogelRefdata() {return  AerogelRefractiveInd;}
  FilterTrData* GetFilterTrData() {return FilterTransData; }
  G4int GetRichTbNumPartEvent() {return  RichTbNumPartEvent; }
  G4int GetRichTbParticleTypeCode() {return RichTbParticleTypeCode; }
  G4int GetRichTbParticleStartPosCode() {return RichTbParticleStartPosCode; } 
  G4int GetRichTbParticleDirectionCode(){return RichTbParticleDirectionCode; } 
  G4int GetRichTbParticleEnergyCode(){return RichTbParticleEnergyCode; } 
  G4double GetRichTbParticleEnergy() {return RichTbParticleEnergy; }
  G4double GetRichTbPhotLowE() {return RichTbPhotLowE; }
  G4double GetRichTbPhotHighE() {return RichTbPhotHighE; }
  G4String GetAerogelRefDataFile()  {return AerogelRefDataInputFileName; }
  G4String GetFilterTransDataFile() { return FilterTransDataInputFileName; }
private:

  G4double CurrentPressureN2;
  G4double CurrentTemperatureN2;
  G4int RefIndexPara;
  G4double MirrorAddTiltX;
  G4double MirrorAddTiltY;
//Graphics setups
  //For Following variables 0 means make the volume invisible;
  //                        1 means make it visible as a solid.
  //                        2 means make it visible as a wireframe.

  G4int RichTbHall_visib; 
  G4int RichTbEnclosure_visib;
  G4int RichTbRadFrame_visib;
  G4int RichTbRadUpW_visib;
  G4int RichTbRadDnW_visib;
  G4int RichTbAerogel_visib;
  G4int RichTbAerogelWrap_visib;
  G4int RichTbFilter_visib;
  G4int RichTbMirror_visib;
  G4int RichTbHpdMaster_visib;
  G4int RichTbHpdEnvelopeTube_visib;
  G4int RichTbHpdQuartzW_visib;
  G4int RichTbHpdPhCathode_visib;
  G4int RichTbHpdSiDet_visib;
  G4int RichTbHpdSectCoat_visib;
  G4int RichTbHpdSiPx_visib;
  G4int RichTbDrawTrajectory_visib;
  G4double HpdPhElectronEnergy;
  G4int FilterTNumber;
  FilterType CurFilterType;
  G4int NumberOfAerogelTiles;
  std::vector<G4int>CurAerogelTNumber;
  std::vector<AerogelType>CurAerogelType;
  G4int WriteOutputFile;
  G4String OutputFileName;
  G4String OutputHistoDirName;
  G4String AerogelRefDataInputFileName;
  G4String FilterTransDataInputFileName;
  G4int WriteMiscOutputFlag;
  G4String MiscOutFileName;
  AerogelRefData* AerogelRefractiveInd;
  FilterTrData* FilterTransData;

  //The following is for particle types.
  // used as the beam.
  //First the number of beam particles per event.
  G4int RichTbNumPartEvent;
  // 0 means piminus, 1 means optical photon,
  // 2 means mixture of piminus and proton.
  G4int RichTbParticleTypeCode;
  //The following poscode = 0 means standard pos (0, 0, -2000.0 mm);
  //                      =1 means special pos implemented.
  G4int  RichTbParticleStartPosCode;
  // The following dire code =0 means dir is 0,0,1. The other
  // codes for getting the beam spread is not implemented yet.
  G4int RichTbParticleDirectionCode;
  G4int RichTbParticleEnergyCode;
  // the following is for charged particle energy in GeV.
 
  G4double RichTbParticleEnergy;
  //The following is for photon energy range in eV for testing optical
  // parameters.
  G4double RichTbPhotLowE;
  G4double RichTbPhotHighE;


};
#endif 






