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
// RichTbMaterialParameters.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "globals.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbMaterialParameters.hh"
#include "FilterTrData.hh"
#include "AerogelTypeSpec.hh"

#include "RichTbAnalysisManager.hh"


void InitializeRichTbMaterial(){

}

std::vector<G4double> InitializePhotonMomentumVector() {

  G4double PhotonEnergyStep=(PhotonMaxEnergy-PhotonMinEnergy)/
                            NumPhotWaveLengthBins;
  std::vector<G4double>PhotMomVect(NumPhotWaveLengthBins);
  for (G4int ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    PhotMomVect[ibin]=PhotonMinEnergy+PhotonEnergyStep*ibin;
  }
  return PhotMomVect;
}
std::vector<G4double> InitN2RefIndex(G4double pressure, G4double temperature){
 
  std::vector<G4double> PmV=InitN2RefPhotW();
  std::vector<G4double> RefN2(NumPhotWaveLengthBins);
  G4double GasRhoN2Cur=GasRhoN2atSTP*(GasTemperature_STP/temperature)*
                                     (pressure/ GasPressure_STP);
  G4double epho,pfe,cpfe;
  for(G4int ibinwn =0; ibinwn<NumPhotWaveLengthBins ; ibinwn++ ){

    epho = PmV[ibinwn]/eV;
    pfe  = SellN2F1/(SellN2E1*SellN2E1 - epho*epho )  +
      SellN2F2/(SellN2E2*SellN2E2 - epho*epho );
    cpfe=0.3738*(GasRhoN2Cur/GasMolWeightN2)*pfe;
    RefN2[ibinwn]=std::pow((1.0+2*cpfe)/(1.0-cpfe),0.5); 
  }
  return RefN2;
}
std::vector<G4double> InitN2RefPhotW() {
  return InitializePhotonMomentumVector() ; 
} 
std::vector<G4double> InitAgelPhotW() {
  return InitializePhotonMomentumVector() ; 
}
std::vector<G4double> InitializeHpdQE(G4int ihpdqe) {
  // Initialize the HPD QE
   G4int iqb;
   if(ihpdqe >= NumHpdTot ) {
     G4cout<<"Wrong HPD Number for QE " <<ihpdqe<<"  vs "
	 <<NumHpdTot <<G4endl;
   }
   std::vector<G4double>qeCurPerCent(NumQEbins);
   if(ihpdqe == 0 ){
    for(iqb=0; iqb<NumQEbins; iqb++){
      qeCurPerCent[iqb] =  Hpd0QEPerCent[iqb]* HpdQEReductionFactor;
   }  
  }
  if(ihpdqe == 1 ){
    for(iqb=0; iqb<NumQEbins; iqb++){
      qeCurPerCent[iqb] =  Hpd1QEPerCent[iqb]* HpdQEReductionFactor;
    }  
  }
  if(ihpdqe == 2 ){
    for(iqb=0; iqb<NumQEbins; iqb++){
      qeCurPerCent[iqb] =  Hpd2QEPerCent[iqb]* HpdQEReductionFactor;
    }  
  }
  if(ihpdqe == 3 ){
    for(iqb=0; iqb<NumQEbins; iqb++){
      qeCurPerCent[iqb] =  Hpd3QEPerCent[iqb]* HpdQEReductionFactor;
    }  
  }

 return  qeCurPerCent;
}
std::vector<G4double> InitializeHpdWaveL(G4int ihpdqe) {
  G4int iqb;
 if(ihpdqe >= NumHpdTot ) {
   G4cout<<"Wrong HPD Number for QE wavelength " <<ihpdqe<<"  vs "
	 <<NumHpdTot <<G4endl;
 }
 // for now all HPDs have the same wavelength bins.
 std::vector<G4double>HpdQEW(NumQEbins);
 for (iqb=0; iqb<NumQEbins; iqb++){
   HpdQEW[iqb]= HpdQEWaveL[iqb];
 }  
 return HpdQEW;
}

void  HistoRichTbMaterialProperties(RichTbRunConfig* RConfig) {


  G4int AerogelNum=0;
  G4double waL=200;

  G4double stepsize=7.0;

  // G4double thickness=(GetCurAerogelLength(AerogelNum))/cm;
  AerogelType CurAerogelType=RConfig-> GetCurAerogelType(AerogelNum);



  G4double Aparam=0.;
  G4double Cparam=0.;
  if(CurAerogelType == AerogelTypeA ) {
    Aparam =   AerogelTypeATotTrans;
    Cparam =  AerogelTypeAClarity*cm/(micrometer*micrometer*micrometer*micrometer);
   }


  
  for(G4int Iabin=0; Iabin<100; Iabin ++ ) {
   
    // G4double waLInmu = waL/1000.0;
    // G4double Aetr = Aparam* std::exp(-Cparam * thickness / std::pow(waLInmu,4) );

    waL += stepsize;
  }






  G4int ihpdqa;
  ihpdqa=0;
  std::vector<G4double>WaveL1 = InitializeHpdWaveL(ihpdqa); 
  std::vector<G4double>QEff1 = InitializeHpdQE(ihpdqa);



}
  

std::vector<G4int> getDeadPixelList(G4int ihpdNum,  G4int){
  std::vector<G4int>DeadPixelList;
  // G4int isc,ipsc;


  if(G4int(DeadPixelList.size()) >  MaxNumDeadPixelPerHpdSect ){
    G4cout<<" Too Many dead Pixels in Hpd "<<DeadPixelList.size()
          <<"   in Hpd "<<ihpdNum<<G4endl;
  }

  return DeadPixelList;
}
std::vector<G4double>GetAerogelRScatLength(AerogelType CurrentAerogelType) {

  std::vector<G4double>AgelRayleighScatLength(NumPhotWaveLengthBins);
  std::vector<G4double>AgelPhotW = InitAgelPhotW();
  G4double aClarity=0.;
  if(CurrentAerogelType == AerogelTypeA ) {
    aClarity=AerogelTypeAClarity/(micrometer*micrometer*micrometer*micrometer);
    
  }else if (CurrentAerogelType == AerogelTypeB ) {
    aClarity=AerogelTypeBClarity/(micrometer*micrometer*micrometer*micrometer);
  }else if (CurrentAerogelType == AerogelTypeC ) {
    aClarity=AerogelTypeCClarity/(micrometer*micrometer*micrometer*micrometer);


  }else if (CurrentAerogelType == AerogelTypeD ) {
    aClarity=AerogelTypeDClarity/(micrometer*micrometer*micrometer*micrometer);

  }else if (CurrentAerogelType == AerogelTypeE ) {
    aClarity=AerogelTypeEClarity/(micrometer*micrometer*micrometer*micrometer);

  }else {G4cout<<"Unknown Aerogel Type for Rayleigh Scat Length "<<G4endl; }

  if(aClarity != 0.0 ) {    
    for(G4int ibinw=0; ibinw<NumPhotWaveLengthBins; ibinw++ ){
      G4double ephoton=AgelPhotW[ibinw]/eV;
      //In the following the 1000 is to convert form nm to micrometer
      G4double wphoton=(PhotMomWaveConv/ephoton)/1000.0;
      AgelRayleighScatLength[ibinw]=(std::pow(wphoton,4))/aClarity;
 
    }
  }

  return  AgelRayleighScatLength;
}
G4double GetCurrentBulkTrans(G4double currentMatRefIndex, 
                          G4double currentNeighbourRefIndex,
                          G4double MaxTotMeasuredTransmission){
  G4double ATrans=MaxTotMeasuredTransmission;
  // G4double ePhot;
  // in the following the energy of the photon is not used since
  // it is only an approximate calulation. 
  G4double na=  currentMatRefIndex;
  G4double nb=  currentNeighbourRefIndex;
  G4double LossAtEntrance=std::pow(((na-nb)/(na+nb)),2.0);
  G4double LossAtExit=std::pow(((nb-na)/(nb+na)),2.0);
  
  G4double LightLossAtExternalSurface= LossAtEntrance+ LossAtExit;
  
  ATrans += LightLossAtExternalSurface;
  if(ATrans >= 1.0) ATrans=1.0;
  return ATrans;
}





