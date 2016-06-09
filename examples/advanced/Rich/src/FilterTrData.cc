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
// Rich advanced example for Geant4
// filterTrData.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <fstream.h>
#include "FilterTrData.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbMaterialParameters.hh" 
FilterTrData::FilterTrData(G4int FilterTypeNum, G4String FilterDataFileName )
   :NumberOfTrBins(NumFilterTransBins),
    TransWaveL(vector<G4double>(NumFilterTransBins)),
    TransValue(vector<G4double>(NumFilterTransBins)),
    TransTotValue(vector<G4double>(NumFilterTransBins))
    {
      // In the G4example the same type of filter is used
      // for all types of filters. 
  FilterTypeNumber= FilterTypeNum;
  FilterTransDataFileName=FilterDataFileName;
  if( FilterTypeNum >= 0 && FilterTypeNum <= 5 ) {
    FilterTypeIndex = GlassD263;  
    FilterRefIndexNominal =  RefIndexGlassD263 ; 
  } else if (FilterTypeNum == -1 ) {
    G4cout<<" No filter  for this run "<<G4endl;
  } else   {G4cout<<"FilterTrData : Unknown Filter Type "<<G4endl; }


  //for now default ref index is provided as that of glassd263.

  CurNeighbourRefIndexNominal= NitrogenNominalRefIndex;
  ReadFilterTrans();
  if(FilterTypeNum >=0 ) {
  FilterThickness = 2.0*FilterHalfZArray[FilterTypeNum];
  }
}
void FilterTrData::ReadFilterTrans() {
  //first get the filter type index
  FilterType Ftype = GetFilterTypeIndex();
  if(Ftype == GlassD263 ) {
      const char* D263FilterTransFile= FilterTransDataFileName.c_str();
      if(D263FilterTransFile != 0 ) {
      G4cout<<"Glass D263 Filter trans is from "<<D263FilterTransFile<<G4endl;
      if(NumberOfTrBins != NumPhotBinGlassD263Trans ) {
        NumberOfTrBins =  NumPhotBinGlassD263Trans;
        TransWaveL.resize(NumberOfTrBins);
        TransValue.resize(NumberOfTrBins);
        TransTotValue.resize(NumberOfTrBins);
      }
      ifstream finpga(D263FilterTransFile);
      G4double wa,trn;
      //set the order to have increasing in wavelength after reading in.
      //the data file is with  decreasing order of wavelength. 
      for (G4int fbin=0; fbin< NumberOfTrBins ; fbin++){
      finpga>>wa;
      finpga>>trn;
      TransWaveL[NumberOfTrBins-fbin-1]=wa;      
//      TransValue[NumberOfTrBins-fbin-1]=trn/100.0;
      TransValue[NumberOfTrBins-fbin-1]=GetCurrentFilTrans(trn);
      TransTotValue[NumberOfTrBins-fbin-1]=GetTotFilTrans(trn);
      }
      }else {
	G4cout<<" Unknown file name for filter transmission Data input "<<G4endl;
        G4cout<<" Please specify the correct file name for filter transmission " <<G4endl;
      }

  } else {
    G4cout<<"Unknown Filter or No filter"<<G4endl;
    G4cout<<" Please provide the Filter transmission Data " <<G4endl;

  }
  

}
G4double FilterTrData::GetCurrentFilTrans(G4double trnData) {
  // No wavelength dependance yet for the Fresnel loss at the
  // surfaces. 
  G4double trn = trnData/100.0;
  G4double na=  FilterRefIndexNominal;
  G4double nb=  CurNeighbourRefIndexNominal;
  G4double LossAtFEntrance=pow(((na-nb)/(na+nb)),2.0);
  G4double LossAtFExit = pow(((nb-na)/(nb+na)),2.0);

  if(trn > 1.0 ) trn=1.0;
  return trn;
}
G4double FilterTrData::GetTotFilTrans(G4double trnData) {
  // No wavelength dependance yet.
 
  return  trnData/100.0; 

}
FilterTrData::~FilterTrData(){ ; }







