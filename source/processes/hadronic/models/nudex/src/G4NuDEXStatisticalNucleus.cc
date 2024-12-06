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
//
// -------------------------------------------------------------------
//
//      Author:        E.Mendoza
// 
//      Creation date: May 2024
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  NuDEX code (https://doi.org/10.1016/j.nima.2022.167894)
// 




#include "G4NuDEXStatisticalNucleus.hh"
#include "G4NuDEXLevelDensity.hh"
#include "G4NuDEXInternalConversion.hh"
#include "G4NuDEXPSF.hh"



G4NuDEXStatisticalNucleus::G4NuDEXStatisticalNucleus(G4int Z,G4int A){

  //The default values for these flags are in "GeneralStatNuclParameters.dat"
  //Can be changed with G4NuDEXStatisticalNucleus::SetSomeInitalParameters(...)
  LevelDensityType=-1;
  PSFflag=-1;
  maxspinx2=-1;
  MinLevelsPerBand=-1;
  BandWidth=0;
  MaxExcEnergy=0;
  BROpt=-1;
  SampleGammaWidths=-1;

  //The default values for these flags are in G4NuDEXStatisticalNucleus::Init(...)
  //Can be changed via G4NuDEXStatisticalNucleus::SetInitialParameters02(...):
  ElectronConversionFlag=-1; // All EC
  KnownLevelsFlag=-1; //Use all known levels
  PrimaryGammasIntensityNormFactor=-1;
  PrimaryGammasEcut=-1; //If primary gammas, do not create new primary gammas going to the "Primary gammas" region.
  Ecrit=-1;
  
  hasBeenInitialized=false;
  NBands=-1;
  theLevels=0;
  theKnownLevels=0;
  NKnownLevels=0; NUnknownLevels=0; NLevels=0; KnownLevelsVectorSize=0;
  theRandom1=0;
  theRandom2=0;
  theRandom3=0;
  theLD=0;
  theICC=0;
  thePSF=0;
  TotalGammaRho=0;
  theThermalCaptureLevelCumulBR=0;
  TotalCumulBR=0;

  Z_Int=Z;
  A_Int=A;

  //Random generators:
  seed1=1234567;
  seed2=1234567;
  seed3=1234567;
  theRandom1= new G4NuDEXRandom(seed1);
  theRandom2= new G4NuDEXRandom(seed2);
  theRandom3= new G4NuDEXRandom(seed3);
  Rand1seedProvided=false; Rand2seedProvided=false; Rand3seedProvided=false;
}

void G4NuDEXStatisticalNucleus::SetInitialParameters02(G4int knownLevelsFlag,G4int electronConversionFlag,G4double primGamNormFactor,G4double primGamEcut,G4double ecrit){
  if(hasBeenInitialized){
    //std::cout<<" ############## Error: G4NuDEXStatisticalNucleus::SetInitialParameters02 cannot be used after initializing the nucleus  ##############"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  if(knownLevelsFlag>=0){KnownLevelsFlag=knownLevelsFlag;}
  if(electronConversionFlag>=0){ElectronConversionFlag=electronConversionFlag;}
  if(primGamNormFactor>=0){PrimaryGammasIntensityNormFactor=primGamNormFactor;}
  if(primGamEcut>=0){PrimaryGammasEcut=primGamEcut;}
  if(ecrit>=0){Ecrit=ecrit;}

}

void G4NuDEXStatisticalNucleus::SetSomeInitalParameters(G4int LDtype,G4int PSFFlag,G4double MaxSpin,G4int minlevelsperband,G4double BandWidth_MeV,G4double maxExcEnergy,G4int BrOption,G4int sampleGammaWidths,unsigned int aseed1,unsigned int aseed2,unsigned int aseed3){

  if(hasBeenInitialized){
    std::cout<<" ############## Error: G4NuDEXStatisticalNucleus::SetSomeInitalParameters cannot be used after initializing the nucleus  ##############"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(LDtype>0){LevelDensityType=LDtype;}
  if(PSFFlag>=0){PSFflag=PSFFlag;}
  if(MaxSpin>0){maxspinx2=(G4int)(2.*MaxSpin+0.01);}
  if(minlevelsperband>0){MinLevelsPerBand=minlevelsperband;}
  if(BandWidth_MeV!=0){BandWidth=BandWidth_MeV;}
  if(maxExcEnergy!=0){MaxExcEnergy=maxExcEnergy;}
  if(BrOption>0){BROpt=BrOption;}
  if(sampleGammaWidths>=0){SampleGammaWidths=sampleGammaWidths;}
  if(aseed1>0){seed1=aseed1; theRandom1->SetSeed(seed1); Rand1seedProvided=true;}
  if(aseed2>0){seed2=aseed2; theRandom2->SetSeed(seed2); Rand2seedProvided=true;}
  if(aseed3>0){seed3=aseed3; theRandom3->SetSeed(seed3); Rand3seedProvided=true;}

}



G4NuDEXStatisticalNucleus::~G4NuDEXStatisticalNucleus(){

  if(theLevels!=0){delete [] theLevels;}
  for(G4int i=0;i<KnownLevelsVectorSize;i++){
    if(theKnownLevels[i].Ndecays>0){
      delete [] theKnownLevels[i].decayFraction;
      delete [] theKnownLevels[i].decayMode;
    }
    if(theKnownLevels[i].NGammas>0){
      delete [] theKnownLevels[i].FinalLevelID;
      delete [] theKnownLevels[i].multipolarity;
      delete [] theKnownLevels[i].Eg;
      delete [] theKnownLevels[i].cumulPtot;
      delete [] theKnownLevels[i].Pg;
      delete [] theKnownLevels[i].Pe;
      delete [] theKnownLevels[i].Icc;
    }
  }
  if(theKnownLevels!=0){delete [] theKnownLevels;}
  if(theRandom1!=0){delete theRandom1;}
  if(theRandom2!=0){delete theRandom2;}
  if(theRandom3!=0){delete theRandom3;}
  if(theLD!=0){delete theLD;}
  if(theICC!=0){delete theICC;}
  if(thePSF!=0){delete thePSF;}
  if(TotalGammaRho!=0){delete [] TotalGammaRho;}
  if(theThermalCaptureLevelCumulBR!=0){delete [] theThermalCaptureLevelCumulBR;}
  if(TotalCumulBR!=0){
    for(G4int i=0;i<NLevels;i++){
      if(TotalCumulBR[i]!=0){delete [] TotalCumulBR[i];}
    }
    delete [] TotalCumulBR;
  }
}


G4int G4NuDEXStatisticalNucleus::Init(const char* dirname,const char* inputfname){

  hasBeenInitialized=true;
  //-------------------------------------------------------------------
  //First, we read data from files:
  G4int check=0;
  char fname[1000],defaultinputfname[1000];
  theLibDir=std::string(dirname);

  //Special (default) input file:
  snprintf(defaultinputfname,1000,"%s/SpecialInputs/ZA_%d.dat",dirname,Z_Int*1000+A_Int);
  G4int HasDefaultInput=ReadSpecialInputFile(defaultinputfname);
  char* definputfn=0;
  if(HasDefaultInput>0){definputfn=defaultinputfname;}
  
  //General statistical parameters:
  snprintf(fname,1000,"%s/GeneralStatNuclParameters.dat",dirname);
  check=ReadGeneralStatNuclParameters(fname); if(check<0){return -1;}

  //Some default, if not initialized yet:
  if(ElectronConversionFlag<0){ElectronConversionFlag=2;} // All EC
  if(KnownLevelsFlag<0){KnownLevelsFlag=1;} //Use all known levels
  if(PrimaryGammasIntensityNormFactor<0){PrimaryGammasIntensityNormFactor=1;}
  if(PrimaryGammasEcut<0){PrimaryGammasEcut=0;} 
  if(Ecrit<0){
    snprintf(fname,1000,"%s/KnownLevels/levels-param.data",dirname);
    check=ReadEcrit(fname); if(check<0){return -1;}
  }

  
  //Level density:
  theLD=new G4NuDEXLevelDensity(Z_Int,A_Int,LevelDensityType);
  check=theLD->ReadLDParameters(dirname,inputfname,definputfn); //if(check<0){return -1;}
  LevelDensityType=theLD->GetLDType(); //because it can be changed by inputfname or due to lack of data
  if(check<0){
    delete theLD; theLD=0;
    Sn=-1; D0=-1; I0=-1000;
  }
  else{
    theLD->GetSnD0I0Vals(Sn,D0,I0);
  }

  //Known level sheme:
  snprintf(fname,1000,"%s/KnownLevels/z%03d.dat",dirname,Z_Int);
  check=ReadKnownLevels(fname); if(check<0){return -1;}  //here we get/crosscheck Sn
  I0=TakeTargetNucleiI0(fname,check); if(check<0){return -1;} //if no I0 --> out

  if(MaxExcEnergy<=0){
    if(Sn>0){
      MaxExcEnergy=Sn-MaxExcEnergy;
    }
    else{
      MaxExcEnergy=1-MaxExcEnergy;
    }
  }

  //If we don't have level density and the known level scheme is not complete, then we can do nothing ...
  if(theLD==0 && Ecrit<MaxExcEnergy){
    std::cout<<" ###### WARNING: No level density and level scheme not complete for ZA="<<1000*Z_Int+A_Int<<" --> Ecrit="<<Ecrit<<" MeV and MaxExcEnergy = "<<MaxExcEnergy<<" MeV ######"<<std::endl;
    return -1;
  }
  //-------------------------------------------------------------------

  //------------------------------------------------------------------- 
  //Init some variables:
  E_unk_min=Ecrit;
  E_unk_max=MaxExcEnergy;
  
  NBands=0;
  if(BandWidth>0){//then we have to create some bands
    NBands=0;
    while(E_unk_min+BandWidth*NBands<MaxExcEnergy){
      NBands++;
    }
    E_unk_max=E_unk_min+BandWidth*NBands;
  }

  Emin_bands=E_unk_min;
  Emax_bands=E_unk_max;
  //-------------------------------------------------------------------


  //Make some checks:
  MakeSomeParameterChecks01();

  //Level scheme:
  //std::cout<<" creating level scheme ..."<<std::endl;
  CreateLevelScheme();
  //std::cout<<" ............. done"<<std::endl;

  if(KnownLevelsFlag==1){
    InsertHighEnergyKnownLevels();
  }

  //We set the precursors for each level:
  for(G4int i=0;i<NLevels;i++){
    theLevels[NLevels-1-i].seed=theRandom2->Integer(4294967295)+1;
  }

  //Internal conversion:
  theICC=new G4NuDEXInternalConversion(Z_Int);
  snprintf(fname,1000,"%s/ICC_factors.dat",dirname);
  theICC->Init(fname);
  theICC->SetRandom4Seed(theRandom3->GetSeed()); //same seed as for generating the cascades

  //PSF:
  thePSF=new G4NuDEXPSF(Z_Int,A_Int);
  thePSF->Init(dirname,theLD,inputfname,definputfn,PSFflag);

  //We compute the missing BR in the known part of the level scheme:
  ComputeKnownLevelsMissingBR();

  //Init TotalGammaRho:
  TotalGammaRho=new G4double[NLevels];
  for(G4int i=0;i<NLevels-1;i++){
    TotalGammaRho[i]=-1;
  }

  //Thermal capture level:
  if(Sn>0 && NLevels>1){
    CreateThermalCaptureLevel();
    GenerateThermalCaptureLevelBR(dirname);
  }

  //Init TotalCumulBR, if BROpt==1,2
  if(BROpt==1 || BROpt==2){
    TotalCumulBR=new G4double*[NLevels];
    for(G4int i=0;i<NLevels;i++){
      TotalCumulBR[i]=0;
    }
  }

  return 0;
}


void G4NuDEXStatisticalNucleus::MakeSomeParameterChecks01(){

  if(LevelDensityType<1 || LevelDensityType>3){
    std::cout<<" ############## Error, LevelDensityType cannot be set to: "<<LevelDensityType<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  
  if(maxspinx2<=0){
    std::cout<<" ############## Error, MaxSpin cannot be set to: "<<maxspinx2/2.<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(MaxExcEnergy<=0){
    std::cout<<" ############## Error, MaxExcEnergy cannot be set to: "<<MaxExcEnergy<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(BROpt<0 || BROpt>2){
    std::cout<<" ############## Error, BROpt cannot be set to: "<<BROpt<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(SampleGammaWidths<0 || SampleGammaWidths>1){
    std::cout<<" ############## Error, SampleGammaWidths cannot be set to: "<<SampleGammaWidths<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  
  if(KnownLevelsFlag!=0 && KnownLevelsFlag!=1){
    std::cout<<" ############## Error, KnownLevelsFlag cannot be set to: "<<KnownLevelsFlag<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  
  if(ElectronConversionFlag!=0 && ElectronConversionFlag!=1 && ElectronConversionFlag!=2){
    std::cout<<" ############## Error, ElectronConversionFlag cannot be set to: "<<ElectronConversionFlag<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(PrimaryGammasIntensityNormFactor<=0){
    std::cout<<" ############## Error, PrimaryGammasIntensityNormFactor cannot be set to: "<<PrimaryGammasIntensityNormFactor<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(PrimaryGammasEcut<0){
    std::cout<<" ############## Error, PrimaryGammasEcut cannot be set to: "<<PrimaryGammasEcut<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(Ecrit<0){
    std::cout<<" ############## Error, Ecrit cannot be set to: "<<Ecrit<<" ##############"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

}


//If InitialLevel==-1 then we start from the thermal capture level
//If ExcitationEnergy>0 then is the excitation energy of the nucleus
//If ExcitationEnergy<0 then is a capture reaction of a neutron with energy -ExcitationEnergy
// return Npar (number of particles emitted). If something goes wrong, returns negative value (for example negative energy transition, which could happen).
G4int G4NuDEXStatisticalNucleus::GenerateCascade(G4int InitialLevel,G4double ExcitationEnergy,std::vector<char>& pType,std::vector<double>& pEnergy,std::vector<double>& pTime){

  pType.clear();
  pEnergy.clear();
  pTime.clear();
  
  if(ExcitationEnergy<0){
    ExcitationEnergy=Sn-(A_Int-1.)/(G4double)A_Int*ExcitationEnergy;
  }
  if(ExcitationEnergy<=0){
    return 0;
  }

  G4int Npar=0;
  G4int f_level=0,multipol=0;
  G4double alpha,E_trans,Exc_ene_i,Exc_ene_f; //icc factor, energy of the transition, initial/final excitation energy
  G4double EmissionTime=0; //in seconds
  G4int NTransition=0;
  //G4double TotalCascadeEnergy1=0,TotalCascadeEnergy2=0;
  
  //Start:
  G4int i_level=InitialLevel;
  Exc_ene_i=ExcitationEnergy;


  if(i_level==0){ //could happen
    pType.push_back('g');
    pEnergy.push_back(Exc_ene_i);
    pTime.push_back(0);
    Npar++;
  }
  
  //Loop:
  while(i_level!=0){

    NTransition++;
    //--------------------------------------------
    //Sample final level:
    if(i_level==-1){ //thermal level
      if(!theThermalCaptureLevelCumulBR){
	f_level=0;
	std::cout<<" ############## NuDEX: WARNING, there are no thermal capture for ZA="<<A_Int+1000*Z_Int-1<<" , with Sn = "<<Sn<<" ##############"<<std::endl; 
      }
      else{
	//Sample final level:
	G4double randnumber=theRandom3->Uniform();
	f_level=-1;
	for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
	  if(theThermalCaptureLevelCumulBR[i]>randnumber){
	    multipol=GetMultipolarity(&theThermalCaptureLevel,&theLevels[i]);
	    f_level=i;
	    break;
	  }
	}
      }
      if(f_level<0){
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      Exc_ene_f=theLevels[f_level].Energy;
    }
    else if(i_level>0){
      f_level=SampleFinalLevel(i_level,multipol,alpha,NTransition);
    }
    else{
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
    //--------------------------------------------

    //Energy of the transition:
    Exc_ene_f=theLevels[f_level].Energy;

    //We sample the final energy if it is a band of levels:
    if(theLevels[f_level].Width!=0){
      Exc_ene_f+=theRandom3->Uniform(-theLevels[f_level].Width,+theLevels[f_level].Width);
    }
    E_trans=Exc_ene_i-Exc_ene_f;
    if(E_trans<=0){
      //std::cout<<"Exc_ene_i = "<<Exc_ene_i<<"  Exc_ene_f = "<<Exc_ene_f<<std::endl;
      //std::cout<<" ####### WARNING: E_trans = "<<E_trans<<" for i="<<i_level<<" with E = "<<theLevels[std::max(i_level,0)].Energy<<" to  f="<<f_level<<" with E = "<<theLevels[f_level].Energy<<" ########"<<std::endl;  
      return -1;
    }
    //------------------------------------------------------------
    //Emission time:
    if(i_level<NKnownLevels && i_level>0){
      if(theKnownLevels[i_level].T12>0){
	EmissionTime+=theRandom3->Exp(theKnownLevels[i_level].T12/std::log(2));
      }
    }
    //------------------------------------------------------------

    //------------------------------------------------------------
    //calculate electron conversion:
    G4bool ele_conv=false;
    if(ElectronConversionFlag>0){
      if(i_level<NKnownLevels && i_level>0){ //ElectronConversionFlag=1,2
	ele_conv=theICC->SampleInternalConversion(E_trans,multipol,alpha); //use the alpha value from the know level value
      }
      else if(ElectronConversionFlag==2){
        ele_conv=theICC->SampleInternalConversion(E_trans,multipol); //calculate alpha (icc factor)
      }
    }
    //------------------------------------------------------------
    //std::cout<<" ---- "<<Exc_ene_i<<"  "<<Exc_ene_f<<"  "<<E_trans<<"  "<<multipol<<"  "<<ele_conv<<std::endl;
    //TotalCascadeEnergy1+=E_trans;
    
    //------------------------------------------------------------
    //Fill result:
    if(ele_conv){
      for(G4int i=0;i<theICC->Ne;i++){
	pType.push_back('e');
	pEnergy.push_back(theICC->Eele[i]);
	pTime.push_back(EmissionTime);
	Npar++;
      }
      for(G4int i=0;i<theICC->Ng;i++){
	pType.push_back('g');
	pEnergy.push_back(theICC->Egam[i]);
	pTime.push_back(EmissionTime);	
	Npar++;
      }
    }
    else{
      pType.push_back('g');
      pEnergy.push_back(E_trans);
      pTime.push_back(EmissionTime);	      
      Npar++;
    }
    //------------------------------------------------------------
    i_level=f_level;
    Exc_ene_i=Exc_ene_f;
  }

  //for(G4int i=0;i<Npar;i++){TotalCascadeEnergy2+=pEnergy[i];}
  //std::cout<<" Total energy: "<<TotalCascadeEnergy1<<"  "<<TotalCascadeEnergy2<<std::endl; getchar();
  
  
  if(i_level!=0){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  return Npar;
}




G4int G4NuDEXStatisticalNucleus::SampleFinalLevel(G4int i_level,G4int& multipolarity,G4double &icc_fac,G4int nTransition){

  if(i_level<=0 || i_level>=NLevels){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  G4double randnumber=theRandom3->Uniform();

  G4int i_knownLevel=-1;
  if(i_level<NKnownLevels){ //then is a known level
    i_knownLevel=i_level;
  }
  if(theLevels[i_level].KnownLevelID>0){ //then is in the unknown part, but we use it as a known level
    if(theKnownLevels[theLevels[i_level].KnownLevelID].NGammas>0){
      i_knownLevel=theLevels[i_level].KnownLevelID;
    }
  }

  if(i_knownLevel>=0){//known part of the spectrum
    theSampledLevel=-1;
    for(G4int j=0;j<theKnownLevels[i_knownLevel].NGammas;j++){
      if(theKnownLevels[i_knownLevel].cumulPtot[j]>randnumber){
	multipolarity=theKnownLevels[i_knownLevel].multipolarity[j];
	icc_fac=theKnownLevels[i_knownLevel].Icc[j];
	return theKnownLevels[i_knownLevel].FinalLevelID[j];
      }
    }
    std::cout<<randnumber<<"  "<<i_knownLevel<<"  "<<theKnownLevels[i_knownLevel].NGammas<<std::endl;
    for(G4int j=0;j<theKnownLevels[i_knownLevel].NGammas;j++){
      std::cout<<theKnownLevels[i_knownLevel].cumulPtot[j]<<std::endl;
    }
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  else{
    icc_fac=-1;
    //------------------------------------------------------------------------------
    //If BROpt==1 or 2, then we store the BR, if not computed, or calculate the final level from it
    if(BROpt==1 || (BROpt==2 && nTransition==1)){
      //maybe the TotalGammaRho[i_level] and BR have not been computed yet:
      if(TotalCumulBR[i_level]==0){
	TotalCumulBR[i_level]=new G4double[i_level];
	TotalGammaRho[i_level]=ComputeDecayIntensities(i_level,TotalCumulBR[i_level]);
      }
      for(G4int j=0;j<i_level;j++){
	if(TotalCumulBR[i_level][j]>randnumber){
	  multipolarity=GetMultipolarity(&theLevels[i_level],&theLevels[j]);
	  return j;
	}
      }
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
    //------------------------------------------------------------------------------

    //BROpt==0
    //------------------------------------------------------------------------------
    // If not, maybe the TotalGammaRho[i_level] has not been computed yet:
    if(TotalGammaRho[i_level]<0){//not computed, we compute it:
      TotalGammaRho[i_level]=ComputeDecayIntensities(i_level);
    }
    theSampledLevel=-1;
    ComputeDecayIntensities(i_level,0,randnumber); // here we compute the final level
    multipolarity=theSampledMultipolarity;
    return theSampledLevel;
    //------------------------------------------------------------------------------
  }

  NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  return 0;
}

void G4NuDEXStatisticalNucleus::ChangeLevelSpinParityAndBR(G4int i_level,G4int newspinx2,G4bool newParity,G4int nlevels,G4double width,unsigned int seed){

  if(i_level==-1){ //change BR of thermal, ignore arguments
    if(Sn>0 && NLevels>1){
      CreateThermalCaptureLevel(seed);
      GenerateThermalCaptureLevelBR(theLibDir.c_str());
    }
    return;
  }

  if(i_level<0 || i_level>=NLevels){
    std::cout<<" i_level = "<<i_level<<" ------ NLevels = "<<NLevels<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  //Do not apply to known levels:
  if(i_level<NKnownLevels || theLevels[i_level].KnownLevelID>0){
    std::cout<<" ####### WARNING: you are trying to change the BR, spin, parity, etc. of a known level --> nothing is done ############"<<std::endl;
    return;
    //NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  theLevels[i_level].spinx2=newspinx2;
  theLevels[i_level].parity=newParity;
  if(seed>0){
    theLevels[i_level].seed=seed;
  }
  else{
    theLevels[i_level].seed=theRandom2->Integer(4294967295)+1;
  }
  if(nlevels>=0){
    theLevels[i_level].NLevels=nlevels;
  }
  if(width>=0){
    theLevels[i_level].Width=width;
  }

  if(TotalGammaRho[i_level]>=0){ //then we have to change TotalGammaRho[i_level]
    G4double* br_vector=0;
    if(TotalCumulBR!=0){
      br_vector=TotalCumulBR[i_level];
    }
    TotalGammaRho[i_level]=ComputeDecayIntensities(i_level,br_vector);
  }

}


//if randnumber<0, return the total TotalGammaRho, and if cumulativeBR!=0, the corresponding cumulativeBR vector is calculated
//if randnumber>0, it is assumed that TotalGammaRho has been already computed 
//          (in the TotalGammaRho[] array or in the TotGR argument) and is used to sample the transition
//     The result is stored in theSampledLevel and theSampledMultipolarity variables
G4double G4NuDEXStatisticalNucleus::ComputeDecayIntensities(G4int i_level,G4double* cumulativeBR,G4double randnumber,G4double TotGR,G4bool AllowE1){

  G4bool  ComputeAlsoBR=false;
  if(cumulativeBR!=0){ComputeAlsoBR=true;}
  //-------------------------------------------------------------------------------------
  //Some checks:
  if(i_level>=NLevels || i_level<0){
    std::cout<<" ############ Error: i = "<<i_level<<" out of range. NLevels = "<<NLevels<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  if(randnumber>0){
    ComputeAlsoBR=false;
    if(TotGR<=0){
      TotGR=TotalGammaRho[i_level];
    }
    if(TotGR<=0){
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
  }
  //-------------------------------------------------------------------------------------

  theRandom2->SetSeed(theLevels[i_level].seed);
  G4double thisTotalGammaRho=0;
  for(G4int j=0;j<i_level;j++){
    //If "solape" then zero:
    if(theLevels[i_level].Energy-theLevels[i_level].Width<=theLevels[j].Energy+theLevels[j].Width){
      thisTotalGammaRho+=0; //not cecessary, but for understanding ...
      if(ComputeAlsoBR){
	cumulativeBR[j]=0;
      }
    }
    else{
      G4double Eg=theLevels[i_level].Energy-theLevels[j].Energy;

      //------------------------------------------------------------------
      //Check which are allowed transitions:
      G4bool E1allowed=true,M1allowed=true,E2allowed=true;
      G4int Lmin=std::abs(theLevels[i_level].spinx2-theLevels[j].spinx2)/2;
      G4int Lmax=(theLevels[i_level].spinx2+theLevels[j].spinx2)/2;
      if(theLevels[i_level].parity==theLevels[j].parity){
	E1allowed=false;
      }
      else{
	M1allowed=false; E2allowed=false;
      }
      if(Lmin>1 || Lmax<1){
	E1allowed=false; M1allowed=false; 
      }
      if(Lmin>2 || Lmax<2){
	E2allowed=false; 
      }
      if(AllowE1){E1allowed=true;}
      //------------------------------------------------------------------

      G4double GammaRho=0,Sumrand2;
      theSampledMultipolarity=-50;
      G4int RealNTransitions=theLevels[i_level].NLevels*theLevels[j].NLevels;

      G4double rand;
      G4double MaxNSamplesForChi2=1000;

      if(E1allowed){
	Sumrand2=RealNTransitions;
	if(SampleGammaWidths==1){ //Porter-Thomas fluctuations
	  Sumrand2=0;
	  if(RealNTransitions>MaxNSamplesForChi2){
	    Sumrand2=RealNTransitions*theRandom2->Gaus(1,std::sqrt(2./RealNTransitions));
	  }
	  else{
	    for(G4int ntr=0;ntr<RealNTransitions;ntr++){
	      rand=theRandom2->Gaus(0,1);
	      Sumrand2+=rand*rand;
	    }
	  }
	}
	GammaRho+=Sumrand2*Eg*Eg*Eg*thePSF->GetE1(Eg,theLevels[i_level].Energy);
	if(randnumber>=0){
	  if(thisTotalGammaRho+GammaRho>=TotGR*randnumber){
	    theSampledMultipolarity=1;
	  }
	}
      }
      if(M1allowed){
	Sumrand2=RealNTransitions;
	if(SampleGammaWidths==1){ //Porter-Thomas fluctuations
	  Sumrand2=0;
	  if(RealNTransitions>MaxNSamplesForChi2){
	    Sumrand2=RealNTransitions*theRandom2->Gaus(1,std::sqrt(2./RealNTransitions));
	  }
	  else{
	    for(G4int ntr=0;ntr<RealNTransitions;ntr++){
	      rand=theRandom2->Gaus(0,1);
	      Sumrand2+=rand*rand;
	    }
	  }
	}
	GammaRho+=Sumrand2*Eg*Eg*Eg*thePSF->GetM1(Eg,theLevels[i_level].Energy);
	if(randnumber>=0){
	  if(thisTotalGammaRho+GammaRho>=TotGR*randnumber){
	    theSampledMultipolarity=-1;
	  }
	}
      }
      if(E2allowed){
	Sumrand2=RealNTransitions;
	if(SampleGammaWidths==1){ //Porter-Thomas fluctuations
	  Sumrand2=0;
	  if(RealNTransitions>MaxNSamplesForChi2){
	    Sumrand2=RealNTransitions*theRandom2->Gaus(1,std::sqrt(2./RealNTransitions));
	  }
	  else{
	    for(G4int ntr=0;ntr<RealNTransitions;ntr++){
	      rand=theRandom2->Gaus(0,1);
	      Sumrand2+=rand*rand;
	    }
	  }
	}
	GammaRho+=Sumrand2*Eg*Eg*Eg*Eg*Eg*thePSF->GetE2(Eg,theLevels[i_level].Energy);
	if(randnumber>=0){
	  if(thisTotalGammaRho+GammaRho>=TotGR*randnumber && theSampledMultipolarity<-10){
	    theSampledMultipolarity=2;
	  }
	}
      }

      /*
      if(i_level==NLevels-1 && j<10){
	std::cout<<j<<"   "<<GammaRho<<"  "<<E1allowed<<"  "<<M1allowed<<"  "<<E2allowed<<"   "<<GammaRho/Eg/Eg/Eg/thePSF->GetE1(Eg,theLevels[i_level].Energy)<<std::endl;
      }
      */

      thisTotalGammaRho+=GammaRho;
      if(ComputeAlsoBR){
	cumulativeBR[j]=GammaRho;
      }
    }

    if(randnumber>=0 && thisTotalGammaRho>=TotGR*randnumber){
      if(theSampledMultipolarity==-50){
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      theSampledLevel=j;
      return -1;
    }
  }

  //If there are no allowed transitions:
  if(randnumber>=0 && thisTotalGammaRho==0){ //if randnumber>0 then TotalGammaRho[i_lev] has been already computed allowing E1 transitions
    return ComputeDecayIntensities(i_level,0,randnumber,TotGR,true);
  }

  if(randnumber>=0){ //then we should not be here
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(thisTotalGammaRho>0){
    if(ComputeAlsoBR){
      for(G4int j=0;j<i_level;j++){
	cumulativeBR[j]/=thisTotalGammaRho;
      }
      for(G4int j=1;j<i_level;j++){
	cumulativeBR[j]+=cumulativeBR[j-1];
      }
      if(std::fabs(cumulativeBR[i_level-1]-1)>1.e-10){
	std::cout<<" ############### Warning:  cumulativeBR["<<i_level<<"]["<<i_level-1<<"]-1 is "<<cumulativeBR[i_level-1]-1<<" ###############"<<std::endl;
      }
    }
  }
  else{
    //std::cout<<" ############### WARNING: total GammaRho for level "<<i_level<<" is "<<thisTotalGammaRho<<". We recalculate it allowing all transitions and assuming the E1 PSF ###############"<<std::endl; 
    //NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    thisTotalGammaRho=ComputeDecayIntensities(i_level,cumulativeBR,-1,-1,true);
  }

  return thisTotalGammaRho;
}


//retrieves the "lowest" allowed multipolarity:
G4int G4NuDEXStatisticalNucleus::GetMultipolarity(Level* theInitialLevel,Level* theFinalLevel){

  if(theInitialLevel->spinx2+theFinalLevel->spinx2==0){
    return 0;
  }
  G4int Lmin=std::abs(theInitialLevel->spinx2-theFinalLevel->spinx2)/2;
  if(Lmin==0){Lmin=1;}
  if(Lmin%2==0){
    if(theInitialLevel->parity==theFinalLevel->parity){
      return +Lmin;
    }
    else{
      return -Lmin;
    }
  }
  else{
    if(theInitialLevel->parity==theFinalLevel->parity){
      return -Lmin;
    }
    else{
      return +Lmin;
    }
  }

  return 0;
}


//Retrieves the closest level of the given spin and parity
//if spinx2<0, then retrieves the closest level
G4int G4NuDEXStatisticalNucleus::GetClosestLevel(G4double Energy,G4int spinx2,G4bool parity){

  //std::cout<<" XXX finding closest level of spin "<<spinx2/2.<<" and parity "<<parity<<" to "<<Energy<<" MeV"<<std::endl;

  //------------------------------------------------------------------------------
  // We try to go closer to the solution, otherwise it takes too much time:
  G4int i_down=0,i_up=NLevels-1;
  G4int i_close_down=0,i_close_up=NLevels-1;
  G4int i_close=0;
  while(i_close_up-i_close_down>10){
    i_close=(i_close_up+i_close_down)/2;
    if(theLevels[i_close].Energy>Energy){
      i_close_up=i_close;
    }
    else{
      i_close_down=i_close;
    }
  }

  for(G4int i=i_close_up;i<NLevels;i++){
    i_up=i;
    if((theLevels[i].spinx2==spinx2 && theLevels[i].parity==parity) || spinx2<0){
      break;
    }
  }
  for(G4int i=i_close_down;i>=0;i--){
    i_down=i;
    if((theLevels[i].spinx2==spinx2 && theLevels[i].parity==parity) || spinx2<0){
      break;
    }
  }
  //------------------------------------------------------------------------------

  G4double MinEnergyDistance=-1,EnergyDistance;
  G4int result=-1;
  for(G4int i=i_down;i<=i_up;i++){
    EnergyDistance=std::fabs(theLevels[i].Energy-Energy);
    if((theLevels[i].spinx2==spinx2 && theLevels[i].parity==parity) || spinx2<0){ //then this is a candidate
      if(EnergyDistance<MinEnergyDistance || MinEnergyDistance<0){
	MinEnergyDistance=EnergyDistance;
	result=i;
      }
    }
  }
  //std::cout<<" XXX found --> "<<result<<std::endl;


  return result;
}


Level* G4NuDEXStatisticalNucleus::GetLevel(G4int i_level){

  if(i_level>=0 && i_level<NLevels){
    return &theLevels[i_level];
  }
  if(i_level==-1){
    return &theThermalCaptureLevel;
  }
  
  std::cout<<" ############ WARNING: for ZA="<<A_Int+1000*Z_Int<<" , requested level i_level="<<i_level<<" does not exist ############"<<std::endl;

  return 0;
}

G4double G4NuDEXStatisticalNucleus::GetLevelEnergy(G4int i_level){

  if(i_level>=0 && i_level<NLevels){
    return theLevels[i_level].Energy;
  }

  return -1;
}


void G4NuDEXStatisticalNucleus::ComputeKnownLevelsMissingBR(){


  for(G4int i=1;i<NKnownLevels;i++){
    if(theKnownLevels[i].NGammas!=0){
      for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
	theKnownLevels[i].multipolarity[j]=GetMultipolarity(&theLevels[i],&theLevels[theKnownLevels[i].FinalLevelID[j]]);
      }
    }
    if(theKnownLevels[i].NGammas==0){
      G4double* tmp=new G4double[i];
      ComputeDecayIntensities(i,tmp);
      for(G4int j=1;j<i;j++){
	if(tmp[j]!=tmp[j-1]){theKnownLevels[i].NGammas++;}
      }
      if(tmp[0]!=0){theKnownLevels[i].NGammas++;}
      if(theKnownLevels[i].NGammas<=0){
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      
      theKnownLevels[i].FinalLevelID=new G4int[theKnownLevels[i].NGammas];
      theKnownLevels[i].multipolarity=new G4int[theKnownLevels[i].NGammas];
      theKnownLevels[i].Eg=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].cumulPtot=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Pg=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Pe=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Icc=new G4double[theKnownLevels[i].NGammas];
      G4int cont=0;
      if(tmp[0]!=0){
	theKnownLevels[i].FinalLevelID[cont]=0;
	theKnownLevels[i].Eg[cont]=theKnownLevels[i].Energy-theKnownLevels[0].Energy;
	theKnownLevels[i].cumulPtot[cont]=tmp[0];
	theKnownLevels[i].multipolarity[cont]=GetMultipolarity(&theLevels[i],&theLevels[theKnownLevels[i].FinalLevelID[cont]]);
	cont++;
      }
      for(G4int j=1;j<i;j++){
	if(tmp[j]!=tmp[j-1]){
	  theKnownLevels[i].FinalLevelID[cont]=j;
	  theKnownLevels[i].Eg[cont]=theKnownLevels[i].Energy-theKnownLevels[j].Energy;
	  theKnownLevels[i].cumulPtot[cont]=tmp[j];
	  theKnownLevels[i].multipolarity[cont]=GetMultipolarity(&theLevels[i],&theLevels[theKnownLevels[i].FinalLevelID[cont]]);
	  cont++;
	}
      }
      delete [] tmp;
      for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
	theKnownLevels[i].Icc[j]=0;
	if(ElectronConversionFlag==2){
	  if(theICC){
	    theKnownLevels[i].Icc[j]=theICC->GetICC(theKnownLevels[i].Eg[j],theKnownLevels[i].multipolarity[j]);
	  }
	}
      }
      G4double alpha=theKnownLevels[i].Icc[0];
      theKnownLevels[i].Pg[0]=theKnownLevels[i].cumulPtot[0]*(1./(alpha+1.));
      theKnownLevels[i].Pe[0]=theKnownLevels[i].cumulPtot[0]*(alpha/(alpha+1.));
      for(G4int j=1;j<theKnownLevels[i].NGammas;j++){
	alpha=theKnownLevels[i].Icc[j];
	theKnownLevels[i].Pg[j]=(theKnownLevels[i].cumulPtot[j]-theKnownLevels[i].cumulPtot[j-1])*(1./(alpha+1.));
	theKnownLevels[i].Pe[j]=(theKnownLevels[i].cumulPtot[j]-theKnownLevels[i].cumulPtot[j-1])*(alpha/(alpha+1.));
      }
    }
  }


}




void G4NuDEXStatisticalNucleus::CreateLevelScheme(){

  //The known levels have been read already
  NLevels=-1;
  Level* theUnkonwnLevels=0;
  if(E_unk_min>=E_unk_max){//Then we know all the level scheme
    NUnknownLevels=0; //will be updated to 1 when creating the capture level
  }
  else{
    //===================================================================
    //Unknown levels:
    G4int maxarraysize=EstimateNumberOfLevelsToFill()*1.1/2.+10000;
    do{
      maxarraysize*=2;
      if(theUnkonwnLevels!=0){delete [] theUnkonwnLevels;}
      //std::cout<<" Max array size of "<<maxarraysize<<std::endl;
      theUnkonwnLevels=new Level[maxarraysize];
      NUnknownLevels=GenerateAllUnknownLevels(theUnkonwnLevels,maxarraysize);
    }while(NUnknownLevels<0);
    //===================================================================
  }

  //===================================================================
  //We define the final level scheme:
  NLevels=NKnownLevels+NUnknownLevels;
  //std::cout<<" There are "<<NLevels<<" levels in total: "<<NKnownLevels<<" known and "<<NUnknownLevels<<" statistically generated "<<std::endl;
  theLevels=new Level[NLevels];
  for(G4int i=0;i<NKnownLevels;i++){
    CopyLevel(&theKnownLevels[i],&theLevels[i]);
  }
  for(G4int i=0;i<NUnknownLevels;i++){
    CopyLevel(&theUnkonwnLevels[i],&theLevels[NKnownLevels+i]);
  }
  delete [] theUnkonwnLevels;
  //===================================================================

  //Final check:
  G4int TotalNIndividualLevels=1;
  for(G4int i=1;i<NLevels;i++){
    TotalNIndividualLevels+=theLevels[i].NLevels;
    if(theLevels[i-1].Energy>theLevels[i].Energy){
      std::cout<<" ########### Error creating the level scheme ###########"<<std::endl;
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
  }

  std::cout<<" NuDEX: Generated statistical nucleus for ZA="<<Z_Int*1000+A_Int<<" up to "<<theLevels[NLevels-1].Energy<<" MeV, with "<<NLevels<<" levels in total: "<<NKnownLevels<<" from the database and "<<NUnknownLevels<<" from statistical models, including bands (without bands --> "<<TotalNIndividualLevels<<" levels). "<<std::endl;
  
}


void G4NuDEXStatisticalNucleus::CreateThermalCaptureLevel(unsigned int seed){

 
  G4int capturespinx2=((std::fabs(I0)+0.5)*2+0.01); //this spin is always possible ...
  G4bool capturepar=true; if(I0<0){capturepar=false;}
  theThermalCaptureLevel.Energy=Sn;
  theThermalCaptureLevel.spinx2=capturespinx2;
  theThermalCaptureLevel.parity=capturepar;
  if(seed>0){
    theThermalCaptureLevel.seed=seed;
  }
  else{
    theThermalCaptureLevel.seed=theRandom2->Integer(4294967295)+1;
  }
  theThermalCaptureLevel.KnownLevelID=-1;
  theThermalCaptureLevel.NLevels=1;
  theThermalCaptureLevel.Width=0;

  NLevelsBelowThermalCaptureLevel=0;
  for(G4int i=0;i<NLevels;i++){
    if(theLevels[i].Energy<theThermalCaptureLevel.Energy){
      NLevelsBelowThermalCaptureLevel++;
    }
  }
  NLevelsBelowThermalCaptureLevel--; // we exclude the last level, transitions to there have no sense

  if(NLevelsBelowThermalCaptureLevel<=0){
    NLevelsBelowThermalCaptureLevel=1; //we can go to the ground state (if Sn>0)
  }
  
}

G4int G4NuDEXStatisticalNucleus::GenerateLevelsInSmallRange(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill){

  //If A_Int even/odd --> spinx2 (spin_val*2) should be even/odd
  if(((A_Int+spinx2)%2)!=0){
    return 0;
  }

  //Get the level density:
  G4double meanNLevels=theLD->Integrate(Emin,Emax,spinx2/2.,parity);

  //Sample total number of levels: ????????
  G4int thisNLevels=0;
  if(meanNLevels>0){
    thisNLevels=(G4int)theRandom1->Poisson(meanNLevels);
  }

  if(thisNLevels>=MaxNLevelsToFill){
    std::cout<<" Warning: not enough space to fill levels "<<std::endl;
    return -1;
  }

  //Distribute the levels in the energy interval: ??????
  for(G4int i=0;i<thisNLevels;i++){
    someLevels[i].Energy=theRandom1->Uniform(Emin,Emax);
    someLevels[i].spinx2=spinx2;
    someLevels[i].parity=parity;
    someLevels[i].seed=0;
    someLevels[i].KnownLevelID=-1;
    someLevels[i].NLevels=1;
    someLevels[i].Width=0;
  }

  return thisNLevels;
}

G4int G4NuDEXStatisticalNucleus::GenerateLevelsInBigRange(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill){

  G4int TotalNLevels=0;
  G4int NIntervals=1000;

  // When the LD changes significantly between two levels the Wigner law does not have sense (ot I don't know how to apply it)
  // So we will sample the levels according to a Poisson Law when rho(E+<D>)/rho(E) is big and according to Wigner when 
  // rho(E+<D>)/rho(E) is small, at higher energies. <D> is the mean energy distance between two levels of the same spin and par.
  // The difference between "big" and "small" is given by:
  G4double WignerRatioThreshold=2;
  G4double LevDenThreshold=1; //if there is less than LevDenThreshold/MeV then Wigner has also no sense

  for(G4int i=0;i<NIntervals;i++){
    G4double emin=Emin+(Emax-Emin)*i/(G4double)NIntervals;
    G4double emax=Emin+(Emax-Emin)*(i+1)/(G4double)NIntervals;
    G4double meanene=(emin+emax)/2.;
    G4double LevDen1=theLD->GetLevelDensity(meanene,spinx2/2.,parity);
    if(LevDen1>LevDenThreshold){
      G4double LevDen2=theLD->GetLevelDensity(meanene+1./LevDen1,spinx2/2.,parity);
      if(LevDen2/LevDen1<WignerRatioThreshold){ //then apply Wigner
	//std::cout<<" Wigner way to generate levels abobe "<<emin<<", being E_unk_min = "<<E_unk_min<<std::endl;
	G4int newExtraLevels=GenerateWignerLevels(emin,Emax,spinx2,parity,&(someLevels[TotalNLevels]),MaxNLevelsToFill-TotalNLevels);
	if(newExtraLevels<0){return -1;}
	TotalNLevels+=newExtraLevels;
	break;
      }
    }
    //then use Poisson:
    G4int newExtraLevels=GenerateLevelsInSmallRange(emin,emax,spinx2,parity,&(someLevels[TotalNLevels]),MaxNLevelsToFill-TotalNLevels);
    if(newExtraLevels<0){return -1;}
    TotalNLevels+=newExtraLevels;
  }

  return TotalNLevels;
}


//Wigner law: p(x)=pi/2*rho*x*exp(-pi/4*rho*rho*x*x), where x is the energy distance between two levels of the same spin and parity
G4int G4NuDEXStatisticalNucleus::GenerateWignerLevels(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill){

  //If A_Int even/odd --> spinx2 (spin_val*2) should be even/odd
  if(((A_Int+spinx2)%2)!=0){
    return 0;
  }

  G4int TotalNLevels=0;

  G4double previousELevel=Emin,nextELevel;
  while(previousELevel<Emax){
    G4double LevDen=theLD->GetLevelDensity(previousELevel,spinx2/2.,parity); //levels/MeV
    G4double arandgamma=theRandom1->Uniform();
    G4double DeltaEMultipliedByLevDen=std::sqrt(-4./3.141592*std::log(1.-arandgamma));
    nextELevel=previousELevel+DeltaEMultipliedByLevDen/LevDen;
    if(nextELevel<Emax){
      someLevels[TotalNLevels].Energy=nextELevel;
      someLevels[TotalNLevels].spinx2=spinx2;
      someLevels[TotalNLevels].parity=parity;
      someLevels[TotalNLevels].seed=0;
      someLevels[TotalNLevels].KnownLevelID=-1;
      someLevels[TotalNLevels].NLevels=1;
      someLevels[TotalNLevels].Width=0;
      TotalNLevels++;
      if(TotalNLevels>=MaxNLevelsToFill){
	std::cout<<" Warning: not enough space to fill levels "<<std::endl;
	//NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	return -1;
      }
    }
    previousELevel=nextELevel;
  }

  return TotalNLevels;

}


//We genereate the levels directly in bands, not individually:
G4int G4NuDEXStatisticalNucleus::GenerateBandLevels(G4int bandmin,G4int bandmax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill){

  //If A_Int even/odd --> spinx2 (spin_val*2) should be even/odd
  if(((A_Int+spinx2)%2)!=0){
    return 0;
  }

  G4double Emin=Emin_bands;
  G4double Emax=Emax_bands;
  G4int TotalNLevels=0;

  if(bandmax>NBands-1){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  for(G4int i=bandmin;i<=bandmax;i++){
    G4double emin=Emin+(Emax-Emin)*i/(G4double)NBands;
    G4double emax=Emin+(Emax-Emin)*(i+1.)/(G4double)NBands;
    G4double AverageNumberOfLevels=theLD->Integrate(emin,emax,spinx2/2.,parity);
    G4int NumberOfLevelsInThisBand=0;
    if(AverageNumberOfLevels>0){
      NumberOfLevelsInThisBand=(G4int)theRandom1->Poisson(AverageNumberOfLevels);
    }
    if(NumberOfLevelsInThisBand>0){
      someLevels[TotalNLevels].Energy=(emax+emin)/2.;
      someLevels[TotalNLevels].spinx2=spinx2;
      someLevels[TotalNLevels].parity=parity;
      someLevels[TotalNLevels].seed=0;
      someLevels[TotalNLevels].KnownLevelID=-1;
      someLevels[TotalNLevels].NLevels=NumberOfLevelsInThisBand;
      someLevels[TotalNLevels].Width=emax-emin;
      TotalNLevels++;
      if(TotalNLevels>=MaxNLevelsToFill){
	std::cout<<" Warning: not enough space to fill levels "<<std::endl;
	//NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	return -1;
      }
    }
  }

  return TotalNLevels;
}



G4int G4NuDEXStatisticalNucleus::GenerateAllUnknownLevels(Level* someLevels,G4int MaxNLevelsToFill){

  G4int TotalNLevels=0,NLev;
  if(E_unk_min>=E_unk_max){return 0;}

  for(G4int spinx2=0;spinx2<=maxspinx2;spinx2++){
    for(G4int ipar=0;ipar<2;ipar++){
      //If A_Int even/odd --> spinx2 (spin_val*2) should be even/odd
      if(((A_Int+spinx2)%2)==0){
	//----------------------------------------------------------------------------------------------------------------------
	G4bool parity=true;
	if(ipar==1){parity=false;}

	//We create random levels between E_unk_min and E_unk_max
	//We will create the levels one by one at low energies and directly in bands at higher energies
	//The limit between the two ranges will be given by E_lim_onebyone
	G4double Emin=E_unk_min;
	G4double Emax=E_unk_max;
	G4double E_lim_onebyone=2.*E_unk_max;
	G4int i_Band_E_lim_onebyone=NBands+1; // band corresponding to E_lim_onebyone

	//----------------------------------------------------
	//Calculate E_lim_onebyone:
#ifndef GENERATEEXPLICITLYALLLEVELSCHEME
	if(NBands>0){
	  if(MinLevelsPerBand<=0){ // All the level scheme in bands
	    E_lim_onebyone=0;
	    i_Band_E_lim_onebyone=0;
	  }
	  else{
	    G4double bandwidth=(Emax_bands-Emin_bands)/NBands;
	    G4double rho_lim_onebyone=3.*(MinLevelsPerBand+10.)/bandwidth; // above this energy we start the creation of bands without sampling the levels one by one
	    E_lim_onebyone=theLD->EstimateInverse(rho_lim_onebyone,spinx2/2.,parity);
	  }
	}
	if(E_unk_max-Emax_bands>0.001){ //then E_unk_max>Emax_bands and we generate all the levels explicitly
	  E_lim_onebyone=2.*E_unk_max;
	  i_Band_E_lim_onebyone=NBands+1;
	}

	// E_lim_onebyone should be in a limit between two bands:
	if(E_lim_onebyone>E_unk_min && E_lim_onebyone<E_unk_max){
	  for(G4int i=0;i<NBands;i++){
	    G4double elow_band=Emin_bands+(Emax_bands-Emin_bands)*i/(G4double)NBands;
	    if(elow_band>E_lim_onebyone){
	      E_lim_onebyone=elow_band;
	      i_Band_E_lim_onebyone=i;
	      break;
	    }
	  }
	}
#endif
	//----------------------------------------------------


	if(E_lim_onebyone>E_unk_min){ //then we have to create some of the levels one by one
	  if(E_lim_onebyone<Emax){
	    Emax=E_lim_onebyone;
	  }
	  NLev=GenerateLevelsInBigRange(Emin,Emax,spinx2,parity,&(someLevels[TotalNLevels]),MaxNLevelsToFill-TotalNLevels);
	  if(NLev<0){return -1;}
	  if(NBands>0 && NLev>0){
	    NLev=CreateBandsFromLevels(NLev,&(someLevels[TotalNLevels]),spinx2,parity);
	  }
	  TotalNLevels+=NLev;
	}

	if(i_Band_E_lim_onebyone<NBands){ //then we have to create some of the levels directly with bands
	  NLev=GenerateBandLevels(i_Band_E_lim_onebyone,NBands-1,spinx2,parity,&(someLevels[TotalNLevels]),MaxNLevelsToFill-TotalNLevels);
	  if(NLev<0){return -1;}
	  TotalNLevels+=NLev;
	}
	//----------------------------------------------------------------------------------------------------------------------
      }
    }
  }


  //Order levels by energy:
  qsort(someLevels,TotalNLevels,sizeof(Level), ComparisonLevels);

  return TotalNLevels;
}

//Junta varios niveles en uno solo, creando distintas bandas, para un spin y paridad determinados.
//Entiende que no hay otros spines ni paridades. Si hay otros, peta.
//Devuelve el numero de niveles actualizado.
G4int G4NuDEXStatisticalNucleus::CreateBandsFromLevels(G4int thisNLevels,Level* someLevels,G4int spinx2,G4bool parity){

  G4double Emin=Emin_bands;
  G4double Emax=Emax_bands;

  Level* theBandLevels=new Level[NBands]; 
  for(G4int i=0;i<NBands;i++){
    G4double emin=Emin+(Emax-Emin)*i/(G4double)NBands;
    G4double emax=Emin+(Emax-Emin)*(i+1.)/(G4double)NBands;
    theBandLevels[i].Energy=(emax+emin)/2.;
    theBandLevels[i].spinx2=spinx2;
    theBandLevels[i].parity=parity;
    theBandLevels[i].seed=0;
    theBandLevels[i].KnownLevelID=-1;
    theBandLevels[i].NLevels=0;
    theBandLevels[i].Width=emax-emin;
    G4int NLevelsInThisBand=0;
    for(G4int j=0;j<thisNLevels;j++){
      if(someLevels[j].spinx2!=spinx2 || someLevels[j].parity!=parity){
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      if(someLevels[j].Energy>=emin && someLevels[j].Energy<=emax){
	NLevelsInThisBand+=someLevels[j].NLevels;
      }
    }
    if(NLevelsInThisBand>=MinLevelsPerBand){
      for(G4int j=0;j<thisNLevels;j++){
	if(someLevels[j].Energy>=emin && someLevels[j].Energy<=emax){
	  theBandLevels[i].NLevels+=someLevels[j].NLevels;
	  someLevels[j].Energy=-1;
	}
      }
    }
  }

  G4int FinalNBands=NBands;

  //Eliminate bands with cero levels:
  for(G4int i=0;i<FinalNBands;i++){
    if(theBandLevels[i].NLevels==0){
      if(i!=FinalNBands-1){
	CopyLevel(&theBandLevels[FinalNBands-1],&theBandLevels[i]);
      }
      i--;
      FinalNBands--;
    }
  }

  //Replace levels with bands and update number of levels:
  G4int NbandsCopied=0;
  for(G4int i=0;i<thisNLevels;i++){
    if(someLevels[i].Energy<0){ 
      if(NbandsCopied<FinalNBands){ //this level is replaced by a band
	CopyLevel(&theBandLevels[NbandsCopied],&someLevels[i]);
	NbandsCopied++;
      }
      else{ //there is no band to replace. Copy the last level here
	CopyLevel(&someLevels[thisNLevels-1],&someLevels[i]);
	i--;
	thisNLevels--;
      }
    }
  }

  //Simple check:
  if(NbandsCopied!=FinalNBands){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  delete [] theBandLevels;
  return thisNLevels;
}


//Esto lo que hace es intentar reemplazar niveles de la parte estadistica por los que se conocen:
G4int G4NuDEXStatisticalNucleus::InsertHighEnergyKnownLevels(){



  G4bool* HasBeenInserted=new G4bool[KnownLevelsVectorSize];
  for(G4int i=0;i<KnownLevelsVectorSize;i++){
    HasBeenInserted[i]=false;
  }

  for(G4int kk=0;kk<2;kk++){ //loop two times: first levels with NGammas>0, then the rest ...
    for(G4int k=1;k<5;k++){ //loop in the distance between levels condition
      G4double MaxEnergyDistance=0.1*k; //The level to replace should be close to it, so first we try with X MeV and then 2X MeV ...
      for(G4int i=NKnownLevels;i<KnownLevelsVectorSize;i++){ //loop in the known levels
	if(theKnownLevels[i].Energy>Sn){break;}
	if(HasBeenInserted[i]==false && (theKnownLevels[i].NGammas>0 || kk>0)){
	  //--------------------------------------------------------------------------------
	  G4double MinEnergyDistance=-1,EnergyDistance=0;
	  G4int thespinx2=theKnownLevels[i].spinx2;
	  G4bool thepar=theKnownLevels[i].parity;
	  G4int unknownLevelID=-1;
	  for(G4int j=NKnownLevels;j<NLevels-1;j++){ // loop in the unknown levels
	    if(theLevels[j].spinx2==thespinx2 && theLevels[j].parity==thepar){
	      EnergyDistance=std::fabs(theKnownLevels[i].Energy-theLevels[j].Energy);
	      if((EnergyDistance<MinEnergyDistance || MinEnergyDistance<0) && EnergyDistance<MaxEnergyDistance && theLevels[j].KnownLevelID<0){
		MinEnergyDistance=EnergyDistance;
		unknownLevelID=j;
	      }
	      //else if(EnergyDistance>MinEnergyDistance && MinEnergyDistance>0){
	      //break;
	      //}
	    }
	  }
	  if(unknownLevelID>0 && theLevels[unknownLevelID].NLevels==1){ //then we replace the stat-level by the known level:
	    //std::cout<<" Copy Level "<<i<<" with E= "<<theKnownLevels[i].Energy<<" into level "<<unknownLevelID<<" with E="<<theLevels[unknownLevelID].Energy<<std::endl;
	    CopyLevel(&theKnownLevels[i],&theLevels[unknownLevelID]);
	    theLevels[unknownLevelID].KnownLevelID=i;
	    HasBeenInserted[i]=true;
	  }
	  //--------------------------------------------------------------------------------
	}
      }
    }
  }
  delete [] HasBeenInserted;

  //We re-order the levels:
  qsort(theLevels,NLevels,sizeof(Level), ComparisonLevels);



  //-----------------------------------------------------------------------------
  //we cannot go to from an inserted level with NGammas>0 to the statistical part of the level scheme (the level ID will then be different in the known and unknown level vectors. There is one exception: if the level has been inserted.
  //so, if it is the case, we change the final level for the closest one:
  for(G4int i=NKnownLevels;i<NLevels;i++){
    if(theLevels[i].KnownLevelID>0){
      G4int knownID=theLevels[i].KnownLevelID;
      for(G4int j=0;j<theKnownLevels[knownID].NGammas;j++){
	if(theKnownLevels[knownID].FinalLevelID[j]>=NKnownLevels){//this cannot be
	  //-----------------------------------------------------
	  G4int i_finalknownlevel=theKnownLevels[knownID].FinalLevelID[j];
	  G4double MinEnergyDistance=-1;
	  G4int i_statlevel=-1;
	  for(G4int k=0;k<i;k++){
	    G4double EnergyDistance=std::fabs(theKnownLevels[i_finalknownlevel].Energy-theLevels[k].Energy);
	    if(EnergyDistance<MinEnergyDistance || MinEnergyDistance<0){
	      MinEnergyDistance=EnergyDistance;
	      i_statlevel=k;
	    }
	  }
	  if(MinEnergyDistance<0){
	    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	  }
	  //std::cout<<" Final known level "<<i_finalknownlevel<<" with E = "<<theKnownLevels[i_finalknownlevel].Energy<<" has been replaced by final level "<<i_statlevel<<" with E = "<<theLevels[i_statlevel].Energy<<std::endl;
	  if(theLevels[i_statlevel].KnownLevelID>=0){ //then is an inserted level, we change only the final level
	    theKnownLevels[knownID].FinalLevelID[j]=i_statlevel;
	  }
	  else{ //is a real statistical level. We change the multipolarity and the alpha:
	    theKnownLevels[knownID].FinalLevelID[j]=i_statlevel;
	    theKnownLevels[knownID].multipolarity[j]=GetMultipolarity(&theLevels[i],&theLevels[i_statlevel]);
	    theKnownLevels[knownID].Eg[j]=theLevels[i].Energy-theLevels[i_statlevel].Energy;
	    theKnownLevels[knownID].Pg[j]=theKnownLevels[knownID].Pg[j]+theKnownLevels[knownID].Pe[j];
	    theKnownLevels[knownID].Pe[j]=0;
	    theKnownLevels[knownID].Icc[j]=0; //set to 0 to avoid problems with the electron conversion
	  }
	  //-----------------------------------------------------
	}
      }
    }
  }
  //-----------------------------------------------------------------------------

  return 0;
}



G4int G4NuDEXStatisticalNucleus::EstimateNumberOfLevelsToFill(){

  if(E_unk_min>=E_unk_max){return 0;}

#ifndef GENERATEEXPLICITLYALLLEVELSCHEME
  if(BandWidth>0){
    return maxspinx2*NBands*MinLevelsPerBand;
  }
#endif

  G4double Emin=E_unk_min;
  G4double Emax=E_unk_max;
  G4double TotalNLevels=0;
  G4double TotalNLevelsInsideBands=0;
  G4double MaxNLevelsWithSameSpinAndParity=0;
  G4double emin,emax,meanEnergy,LevDen,meanNLevels;
  G4int NIntervals=1000;
  for(G4int spinx2=0;spinx2<=maxspinx2;spinx2++){
    G4double TotalNLevesWithSameSpin=0;
    for(G4int i=0;i<NIntervals;i++){
      emin=Emin+(Emax-Emin)*i/(G4double)NIntervals;
      emax=Emin+(Emax-Emin)*(i+1)/(G4double)NIntervals;
      meanEnergy=(emax+emin)/2.;

      //Positive parity:
      LevDen=theLD->GetLevelDensity(meanEnergy,spinx2/2.,true); //levels/MeV
      meanNLevels=LevDen*(emax-emin); 
      TotalNLevels+=meanNLevels;
      TotalNLevesWithSameSpin+=meanNLevels;
      if(NBands>0 && meanEnergy>Emin_bands && meanEnergy<Emax_bands){ //if there are bands and these levels go inside:
	TotalNLevelsInsideBands+=meanNLevels;
      }


      //Negative parity:
      LevDen=theLD->GetLevelDensity(meanEnergy,spinx2/2.,false); //levels/MeV
      meanNLevels=LevDen*(emax-emin);
      TotalNLevels+=meanNLevels;
      TotalNLevesWithSameSpin+=meanNLevels;
      if(NBands>0 && meanEnergy>Emin_bands && meanEnergy<Emax_bands){ //if there are bands and these levels go inside:
	TotalNLevelsInsideBands+=meanNLevels;
      }

    }
    if(MaxNLevelsWithSameSpinAndParity<TotalNLevesWithSameSpin){MaxNLevelsWithSameSpinAndParity=TotalNLevesWithSameSpin;}
  }

  MaxNLevelsWithSameSpinAndParity/=2.;

  if(TotalNLevelsInsideBands>0){
    //then there are some levels inside bands and the amount of levels needed will be reduced:
    G4double TotalNLevelsOutsideBands=TotalNLevels-TotalNLevelsInsideBands;
    G4double MaxNLevelsNeededForBands=NBands*maxspinx2*MinLevelsPerBand;
    TotalNLevels=MaxNLevelsWithSameSpinAndParity+TotalNLevelsOutsideBands+MaxNLevelsNeededForBands;
  }

  return (G4int)TotalNLevels;
}


G4double G4NuDEXStatisticalNucleus::TakeTargetNucleiI0(const char* fname,G4int& check){

  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  check=0;
  
  char buffer[1000];
  G4int aZ,aA;
  while(in.get(buffer,6)){
    in.get(buffer,6); aA=atoi(buffer); 
    in.get(buffer,6); aZ=atoi(buffer);
    if(aZ==Z_Int && aA==A_Int-1){
      break;
    }
    in.ignore(10000,'\n');
  }
  if(!in.good()){
    in.close();
    check=-1;
    //NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  in.ignore(10000,'\n');
  G4double spin,par;
  in.get(buffer,16);
  in.get(buffer,6);   spin=std::fabs(atof(buffer)); // some spins are negative ???
  in.get(buffer,4);   par=atof(buffer); 

  in.close();

  if(par<0){return -spin;}

  return spin;
}

G4double G4NuDEXStatisticalNucleus::ReadKnownLevels(const char* fname){


  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  char buffer[1000];
  G4int aZ,aA;
  while(in.get(buffer,6)){
    in.get(buffer,6); aA=atoi(buffer); 
    in.get(buffer,6); aZ=atoi(buffer); 
    if(aZ==Z_Int && aA==A_Int){
      in.get(buffer,6); KnownLevelsVectorSize=atoi(buffer);
      in.get(buffer,16);
      in.get(buffer,13); G4double newSn=atof(buffer);
      if(Sn>0 && std::fabs(Sn-newSn)>0.01){
	std::cout<<" ######## WARNING: Sn value from the level density file ("<<Sn<<") is different than the one from the known levels file ("<<newSn<<"). We use the first value. ########"<<std::endl;
      }
      else if(Sn<0){
	Sn=newSn;
      }
      break;
    }
    in.ignore(10000,'\n');
  }

  if(!in.good()){
    in.close(); return -1;
  }

  in.ignore(10000,'\n');

  NKnownLevels=0;
  theKnownLevels=new KnownLevel[KnownLevelsVectorSize];
  for(G4int i=0;i<KnownLevelsVectorSize;i++){theKnownLevels[i].NGammas=0;}
  G4double spin,par;
  for(G4int i=0;i<KnownLevelsVectorSize;i++){
    in.get(buffer,4);   theKnownLevels[i].id=atoi(buffer)-1; 
    in.get(buffer,2); 
    in.get(buffer,11);  theKnownLevels[i].Energy=atof(buffer); 
    in.get(buffer,2); 
    in.get(buffer,6);   spin=atof(buffer); 
    in.get(buffer,4);   par=atof(buffer); 
    if((spin<0 || par==0) && theKnownLevels[i].Energy<Ecrit){
      std::cout<<" ######## WARNING: Spin and parity for level "<<i<<" is s="<<spin<<" p="<<par<<" for Z="<<Z_Int<<", A="<<A_Int<<" ########"<<std::endl;
      if(spin<0){
	spin=0;
	if(i>1){ //Random spin, same as one of the levels below this one:
	  G4int i_sampleSpin=theRandom1->Integer(i-1);
	  spin=theKnownLevels[i_sampleSpin].spinx2/2.;
	}
      }
      if(par==0){
	par=1;
	if(theRandom1->Uniform(-1,1)<0){par=-1;}
      }
    }
    in.get(buffer,2); 
    in.get(buffer,11);  theKnownLevels[i].T12=atof(buffer); 
    in.get(buffer,4);   theKnownLevels[i].NGammas=atoi(buffer);
    if(theKnownLevels[i].NGammas>0){
      if(spin<0){
	spin=0;
	if(i>1){ //Random spin, same as one of the levels below this one:
	  G4int i_sampleSpin=theRandom1->Integer(i-1);
	  spin=theKnownLevels[i_sampleSpin].spinx2/2.;
	}
      }
      if(par==0){
	par=1;
	if(theRandom1->Uniform(-1,1)<0){par=-1;}
      }
    }    
    theKnownLevels[i].spinx2=(G4int)(spin*2+0.01);
    if(par>0){theKnownLevels[i].parity=true;}else{theKnownLevels[i].parity=false;}

    //---------------------------------
    //decay modes:
    in.get(buffer,27); //dummy
    in.get(buffer,4);
    G4int decays=theKnownLevels[i].Ndecays=atoi(buffer);
    if(decays>0){
      theKnownLevels[i].decayFraction=new G4double[decays];
      theKnownLevels[i].decayMode=new std::string[decays];
    }
    for(G4int j=0;j<decays;j++){
      in.get(buffer,5);
      in.get(buffer,11);
      theKnownLevels[i].decayFraction[j]=atof(buffer);
      in.get(buffer,2);
      in.get(buffer,8);
      theKnownLevels[i].decayMode[j]=std::string(buffer);
    }
    //----------------------------------

    in.ignore(10000,'\n');

    if(theKnownLevels[i].NGammas>0){
      theKnownLevels[i].FinalLevelID=new G4int[theKnownLevels[i].NGammas];
      theKnownLevels[i].multipolarity=new G4int[theKnownLevels[i].NGammas];
      theKnownLevels[i].Eg=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].cumulPtot=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Pg=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Pe=new G4double[theKnownLevels[i].NGammas];
      theKnownLevels[i].Icc=new G4double[theKnownLevels[i].NGammas];
    }

    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      in.get(buffer,40); 

      in.get(buffer,5);  theKnownLevels[i].FinalLevelID[j]=atoi(buffer)-1; 
      theKnownLevels[i].multipolarity[j]=0;
      in.get(buffer,2); 
      in.get(buffer,11);  theKnownLevels[i].Eg[j]=atof(buffer); 
      in.get(buffer,2); 
      in.get(buffer,11);  theKnownLevels[i].Pg[j]=atof(buffer); 
      in.get(buffer,2); 
      in.get(buffer,11);  theKnownLevels[i].Pe[j]=atof(buffer); 
      in.get(buffer,2); 
      in.get(buffer,11);  theKnownLevels[i].Icc[j]=atof(buffer); 
      theKnownLevels[i].cumulPtot[j]=theKnownLevels[i].Pg[j]*(1+theKnownLevels[i].Icc[j]); //we rely in Pg and Icc, where Icc=Pe/Pg
      in.ignore(10000,'\n');
      if(theKnownLevels[i].FinalLevelID[j]>=i+1){
	std::cout<<" ######## Error reading file "<<fname<<" ########"<<std::endl;
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
    }

    if(theKnownLevels[i].id!=i || !in.good()){
      std::cout<<" ######## Error reading file "<<fname<<" ########"<<std::endl;
      std::cout<<" Level "<<i<<" has id = "<<theKnownLevels[i].id<<std::endl;
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }

    //--------------------------------------------------------------------------------------
    //normalize, and put cumulPtot as cumulative. Re-calculate Pe
    G4double totalEmissionProb=0;
    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      totalEmissionProb+=theKnownLevels[i].cumulPtot[j];
    }
    //------------------------------------
    if(totalEmissionProb==0){//sometimes all the levels have 0 prob.
      for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
	theKnownLevels[i].cumulPtot[j]=1./theKnownLevels[i].NGammas;
	theKnownLevels[i].Pg[j]=theKnownLevels[i].cumulPtot[j]/(1.+theKnownLevels[i].Icc[j]);
      }
      totalEmissionProb=1;
    }
    //------------------------------------
    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      theKnownLevels[i].cumulPtot[j]/=totalEmissionProb;
      theKnownLevels[i].Pg[j]/=totalEmissionProb;
      theKnownLevels[i].Pe[j]=theKnownLevels[i].Icc[j]*theKnownLevels[i].Pg[j];
    }
    for(G4int j=1;j<theKnownLevels[i].NGammas;j++){
      theKnownLevels[i].cumulPtot[j]+=theKnownLevels[i].cumulPtot[j-1];
    }
    //--------------------------------------------------------------------------------------

    if(theKnownLevels[i].Energy<=Ecrit*1.0001){
      NKnownLevels++;
    }
    /*
    else{
      break;
    }
    */
  }

  if(!in.good()){
    in.close(); return -1;
  }

  in.close();

  return 0;
}



G4double G4NuDEXStatisticalNucleus::ReadEcrit(const char* fname){

  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  G4int aZ,aA;
  char word[200];
  Ecrit=-1;
  for(G4int i=0;i<4;i++){in.ignore(10000,'\n');}
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      in>>word>>word>>word>>word>>word>>word>>word>>word>>word>>Ecrit; break;
    }
    in.ignore(10000,'\n');
  }
  in.close();

  return Ecrit;
}

G4int G4NuDEXStatisticalNucleus::ReadSpecialInputFile(const char* fname){

  std::string word;
  std::ifstream in(fname);
  if(!in.good()){
    in.close();
    return -1;
  }
  G4double MaxSpin;
  while(in>>word){
    if(word.c_str()[0]=='#'){in.ignore(10000,'\n');}
    if(word==std::string("END")){break;}
    //now we will take values only if they have not been set yet:
    else if(word==std::string("LEVELDENSITYTYPE")){if(LevelDensityType<0){in>>LevelDensityType;}}
    else if(word==std::string("MAXSPIN")){if(maxspinx2<0){in>>MaxSpin; maxspinx2=(G4int)(2.*MaxSpin+0.01);}}
    else if(word==std::string("MINLEVELSPERBAND")){if(MinLevelsPerBand<0){in>>MinLevelsPerBand;}}
    else if(word==std::string("BANDWIDTH_MEV")){if(BandWidth==0){in>>BandWidth;}}
    else if(word==std::string("MAXEXCENERGY_MEV")){if(MaxExcEnergy==0){in>>MaxExcEnergy;}}
    else if(word==std::string("ECRIT_MEV")){if(Ecrit<0){in>>Ecrit;}}
    else if(word==std::string("KNOWNLEVELSFLAG")){if(KnownLevelsFlag<0){in>>KnownLevelsFlag;}}

    else if(word==std::string("PSF_FLAG")){if(PSFflag<0){in>>PSFflag;}}
    else if(word==std::string("BROPTION")){if(BROpt<0){in>>BROpt;}}
    else if(word==std::string("SAMPLEGAMMAWIDTHS")){if(SampleGammaWidths<0){in>>SampleGammaWidths;}}
      
    else if(word==std::string("ELECTRONCONVERSIONFLAG")){if(ElectronConversionFlag<0){in>>ElectronConversionFlag;}}
    else if(word==std::string("PRIMARYTHCAPGAMNORM")){if(PrimaryGammasIntensityNormFactor<0){in>>PrimaryGammasIntensityNormFactor;}}
    else if(word==std::string("PRIMARYGAMMASECUT")){if(PrimaryGammasEcut<0){in>>PrimaryGammasEcut;}}
  }
  in.close();
  return 1;
}

G4int G4NuDEXStatisticalNucleus::ReadGeneralStatNuclParameters(const char* fname){

  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  char line[1000];
  in.getline(line,1000);
  in.getline(line,1000);
  G4int tmpZ,tmpA,tmpLDtype,tmpPSFflag,tmpMaxSpin,tmpMinLev,tmpBrOption,tmpSampleGW;
  G4double tmpBandWidth,tmpMaxExcEnergy;
  unsigned int tmpseed1,tmpseed2,tmpseed3;
  G4int finalLDtype=0,finalPSFflag=0,finalMaxSpin=0,finalMinLev=0,finalBrOption=0,finalSampleGW=0;
  G4double finalBandWidth=0,finalMaxExcEnergy=0;
  unsigned int finalseed1=0,finalseed2=0,finalseed3=0;
  G4bool NuclDataHasBeenRead=false;
  G4bool DefaulDataHasBeenRead=false;
  while(in>>tmpZ>>tmpA>>tmpLDtype>>tmpPSFflag>>tmpMaxSpin>>tmpMinLev>>tmpBandWidth>>tmpMaxExcEnergy>>tmpBrOption>>tmpSampleGW>>tmpseed1>>tmpseed2>>tmpseed3){
    G4bool TakeData=false;
    if(tmpZ==Z_Int && tmpA==A_Int){ //then this is our nucleus
      NuclDataHasBeenRead=true;
      TakeData=true;
    }
    else if(tmpZ==0 && tmpA==0 && NuclDataHasBeenRead==false){ //default, only if our nucleus has not been read
      DefaulDataHasBeenRead=true;
      TakeData=true;
    }
    if(TakeData){ 
      finalLDtype=tmpLDtype; finalPSFflag=tmpPSFflag; finalMaxSpin=tmpMaxSpin; finalMinLev=tmpMinLev; finalBrOption=tmpBrOption; finalSampleGW=tmpSampleGW;
      finalBandWidth=tmpBandWidth; finalMaxExcEnergy=tmpMaxExcEnergy;
      finalseed1=tmpseed1; finalseed2=tmpseed2; finalseed3=tmpseed3;
    }
  }
  in.close();

  if(NuclDataHasBeenRead==false && DefaulDataHasBeenRead==false){
    std::cout<<" ######## Error reading "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  // Replace data only if it has not been set by the SetSomeInitalParameters method:
  if(LevelDensityType<0){LevelDensityType=finalLDtype;}
  if(PSFflag<0){PSFflag=finalPSFflag;}
  if(maxspinx2<0){maxspinx2=(G4int)(2.*finalMaxSpin+0.01);}
  if(MinLevelsPerBand<0){MinLevelsPerBand=finalMinLev;}
  if(BandWidth==0){BandWidth=finalBandWidth;}
  if(MaxExcEnergy==0){MaxExcEnergy=finalMaxExcEnergy;}
  if(BROpt<0){BROpt=finalBrOption;}
  if(SampleGammaWidths<0){SampleGammaWidths=finalSampleGW;}
  if(Rand1seedProvided==false){seed1=finalseed1; theRandom1->SetSeed(finalseed1);}
  if(Rand2seedProvided==false){seed2=finalseed2; theRandom2->SetSeed(finalseed2);}
  if(Rand3seedProvided==false){seed3=finalseed3; theRandom3->SetSeed(finalseed3);}

  //-------------------------------------------------------
  //Now some checks:
  if(maxspinx2<1){
    std::cout<<" ######## Error: maximum spin for generating the statistical nucleus with A="<<A_Int<<" and Z="<<Z_Int<<" has been set to "<<maxspinx2/2.<<" ########"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  if(MinLevelsPerBand<=0 && BandWidth<0){
    std::cout<<" ######## Error: MinLevelsPerBand and BandWidth for generating the statistical nucleus with A="<<A_Int<<" and Z="<<Z_Int<<" has been set to MinLevelsPerBand="<<MinLevelsPerBand<<" and BandWidth="<<BandWidth<<" ########"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  if(BROpt!=0 && BROpt!=1 && BROpt!=2){
    std::cout<<" ######## Error: BROpt for generating the statistical nucleus with A="<<A_Int<<" and Z="<<Z_Int<<" has been set to BROpt="<<BROpt<<", and has to be BROpt=0,1 or 2 ########"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  if(SampleGammaWidths!=0 && SampleGammaWidths!=1){
    std::cout<<" ######## Error: SampleGammaWidths parameter for generating the statistical nucleus with A="<<A_Int<<" and Z="<<Z_Int<<" has been set to SampleGammaWidths="<<SampleGammaWidths<<", and has to be SampleGammaWidths=0 or 1 ########"<<std::endl; NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  //-------------------------------------------------------


  return 0;
}


void G4NuDEXStatisticalNucleus::GenerateThermalCaptureLevelBR(const char* dirname){

  char fname[1000];
  snprintf(fname,1000,"%s/PrimaryCaptureGammas.dat",dirname);

  G4int aZA,ng=0;
  char word[200];
  G4double ELevel;
  G4double ThEg[1000],ThI[1000];
  //G4double TotalThI=0;

  //We read the gamma intensities from the file, if they are there:
  std::ifstream in(fname);
  while(in>>word){
    if(word[0]=='Z' && word[1]=='A'){
      in>>aZA;
      if(aZA==Z_Int*1000+A_Int){
	in.ignore(10000,'\n');
	in>>ELevel>>ng;
	in.ignore(10000,'\n');
	ELevel*=1.e-3; // keV to MeV
	//Check if ELevel is very close to Sn:
	if(ELevel/Sn>1.001 || ELevel/Sn<0.999){
	  std::cout<<" ########## WARNING: ELevel = "<<ELevel<<" and Sn = "<<Sn<<" for ZA = "<<aZA<<" ##########"<<std::endl;
	}
	for(G4int i=0;i<ng;i++){
	  in>>ThEg[i]>>ThI[i];
	  ThEg[i]/=1.e3; // keV to MeV
	  ThI[i]/=100.;  // percent to no-percent
	  ThI[i]*=PrimaryGammasIntensityNormFactor; //renormalization
	  //TotalThI+=ThI[i];
	}
	break;
      }
    }
    in.ignore(10000,'\n');
  }
  in.close();

  if(theThermalCaptureLevelCumulBR){delete [] theThermalCaptureLevelCumulBR;}
  theThermalCaptureLevelCumulBR=new G4double[NLevelsBelowThermalCaptureLevel]; //the final result
  for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
    theThermalCaptureLevelCumulBR[i]=0;
  }

  //------------------------------------------------------------------------------------------------------
  //We calculate which transitrions go to an existing level and the total intensity of the transitions:
  G4double totalThGInt=0,ENDSFLevelEnergy=0,MinlevelDist=0,MinlevelDist_known=0,LevDist=0;
  G4int i_closest=0,i_closest_known=0;
  G4double MaxAllowedLevelDistance=0.010; //10 keV
  G4bool ComputePrimaryGammasEcut=false;
  if(PrimaryGammasEcut==0){ComputePrimaryGammasEcut=true;}
  
  //We take only those intensities going to our levels:
  for(G4int i=0;i<ng;i++){
    ENDSFLevelEnergy=ELevel-ThEg[i];
    if(ComputePrimaryGammasEcut && PrimaryGammasEcut<ENDSFLevelEnergy){
      PrimaryGammasEcut=ENDSFLevelEnergy;
    }
    i_closest=0;
    MinlevelDist=1.e20;
    i_closest_known=0;
    MinlevelDist_known=1.e20;
    for(G4int j=0;j<NLevelsBelowThermalCaptureLevel;j++){
      LevDist=std::fabs(ENDSFLevelEnergy-theLevels[j].Energy);
      if(theLevels[j].KnownLevelID>=0){ //then this is a known level. We priorize known levels.
	if(LevDist<MinlevelDist_known){
	  MinlevelDist_known=LevDist;
	  i_closest_known=j;
	}
      }
      if(LevDist<MinlevelDist){
	MinlevelDist=LevDist;
	i_closest=j;
      }
    }
    if(MinlevelDist_known<MaxAllowedLevelDistance){ // We priorize known levels.
      theThermalCaptureLevelCumulBR[i_closest_known]=ThI[i];
      totalThGInt+=theThermalCaptureLevelCumulBR[i_closest_known];
    }
    else if(MinlevelDist<MaxAllowedLevelDistance){
      theThermalCaptureLevelCumulBR[i_closest]=ThI[i];
      totalThGInt+=theThermalCaptureLevelCumulBR[i_closest];
    }
  }
  //if(TotalThI>0){std::cout<<" NuDEX: Primary thermal gammas for ZA="<<Z_Int*1000+A_Int<<" file:  "<<TotalThI<<" accepted:  "<<totalThGInt<<" ratio: "<<totalThGInt/TotalThI*100.<<" %"<<std::endl;}
  //else{std::cout<<" Primary thermal gammas for "<<Z_Int*1000+A_Int<<" file:  "<<TotalThI<<std::endl;}
  std::cout<<" NuDEX: Primary thermal gammas for ZA="<<Z_Int*1000+A_Int<<" found in the database: "<<totalThGInt*100.<<" %"<<std::endl;
  //------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------------
  //If we don't have all the info, we take the rest of the BR from one of the levels, but with a proper normalization:
  if(totalThGInt<0.95){
    G4double TotalNeededIntensity=1.-totalThGInt;
    G4double oldInt,TotalOldIntensity=0;
    //-------------------------
    //We put the capture level replacing one of the levels and calculate the BR:
    Level tmpLevel;
    CopyLevel(&theLevels[NLevelsBelowThermalCaptureLevel],&tmpLevel);
    CopyLevel(&theThermalCaptureLevel,&theLevels[NLevelsBelowThermalCaptureLevel]);
    G4double* CumulBR_th_v2=new G4double[NLevelsBelowThermalCaptureLevel];
    ComputeDecayIntensities(NLevelsBelowThermalCaptureLevel,CumulBR_th_v2);
    CopyLevel(&tmpLevel,&theLevels[NLevelsBelowThermalCaptureLevel]);
    //-------------------------
    
    for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
      if(i==0){oldInt=CumulBR_th_v2[i];}else{oldInt=CumulBR_th_v2[i]-CumulBR_th_v2[i-1];}
      if(theThermalCaptureLevelCumulBR[i]==0 && theLevels[i].Energy>=PrimaryGammasEcut){TotalOldIntensity+=oldInt;}
    }

    if(TotalOldIntensity>0){
      for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
	if(theThermalCaptureLevelCumulBR[i]==0 && theLevels[i].Energy>=PrimaryGammasEcut){
	  if(i==0){oldInt=CumulBR_th_v2[i];}else{oldInt=CumulBR_th_v2[i]-CumulBR_th_v2[i-1];}
	  theThermalCaptureLevelCumulBR[i]=oldInt*TotalNeededIntensity/TotalOldIntensity;
	}
      }
    }
    delete [] CumulBR_th_v2;
  }
  //------------------------------------------------------------------------------------------------------

  //theThermalCaptureLevelCumulBR still not cumulative

  for(G4int i=1;i<NLevelsBelowThermalCaptureLevel;i++){
    theThermalCaptureLevelCumulBR[i]+=theThermalCaptureLevelCumulBR[i-1];
  }
  for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
   theThermalCaptureLevelCumulBR[i]/=theThermalCaptureLevelCumulBR[NLevelsBelowThermalCaptureLevel-1];
  }

}

//Cambiamos las intensidades de los "primary gammas". Del correspondiente al ninvel con energía "LevelEnergy"
void G4NuDEXStatisticalNucleus::ChangeThermalCaptureLevelBR(G4double LevelEnergy,G4double absoluteIntensity){

  if(!theThermalCaptureLevelCumulBR){return;}
  G4int level_id=GetClosestLevel(LevelEnergy,-1,true);
  if(level_id<0 || level_id>=NLevelsBelowThermalCaptureLevel){
    std::cout<<" ############## WARNING in "<<__FILE__<<", line "<<__LINE__<<" ##############"<<std::endl;
    std::cout<<"  ---> "<<level_id<<"  "<<LevelEnergy<<std::endl;
  }

  for(G4int i=NLevelsBelowThermalCaptureLevel-1;i>0;i--){
    theThermalCaptureLevelCumulBR[i]-=theThermalCaptureLevelCumulBR[i-1];
  }
  G4double OldIntensity=theThermalCaptureLevelCumulBR[level_id];
  theThermalCaptureLevelCumulBR[level_id]=absoluteIntensity*(1.-OldIntensity)/(1.-absoluteIntensity);

  for(G4int i=1;i<NLevelsBelowThermalCaptureLevel;i++){
    theThermalCaptureLevelCumulBR[i]+=theThermalCaptureLevelCumulBR[i-1];
  }
  for(G4int i=0;i<NLevelsBelowThermalCaptureLevel;i++){
   theThermalCaptureLevelCumulBR[i]/=theThermalCaptureLevelCumulBR[NLevelsBelowThermalCaptureLevel-1];
  }
  if(level_id==0){
    std::cout<<" Thermal primary gammas to level "<<level_id<<", with E="<<theLevels[level_id].Energy<<" MeV changed from "<<OldIntensity<<" to "<<theThermalCaptureLevelCumulBR[level_id]<<std::endl;
  }
  else{
    std::cout<<" Thermal primary gammas to level "<<level_id<<", with E="<<theLevels[level_id].Energy<<" MeV changed from "<<OldIntensity<<" to "<<theThermalCaptureLevelCumulBR[level_id]-theThermalCaptureLevelCumulBR[level_id-1]<<std::endl;
  }
}

void G4NuDEXStatisticalNucleus::PrintParameters(std::ostream &out){

  out<<" ###################################################################################### "<<std::endl;
  out<<" GENERAL_PARS"<<std::endl;
  out<<" Z = "<<Z_Int<<"  A = "<<A_Int<<std::endl;
  out<<" Sn = "<<Sn<<"  I0(ZA-1) = "<<I0<<std::endl;
  if(theLD!=0){theLD->PrintParameters(out);}
  else{out<<" No level density"<<std::endl;}
  out<<" PSFflag = "<<PSFflag<<std::endl;
  out<<" Ecrit = "<<Ecrit<<std::endl;
  out<<" E_unknown_min = "<<E_unk_min<<"  E_unknown_max = "<<E_unk_max<<std::endl;
  out<<" maxspin = "<<maxspinx2/2.<<std::endl;
  out<<" MaxExcEnergy = "<<MaxExcEnergy<<std::endl;
  out<<" NBands = "<<NBands<<"  MinLevelsPerBand = "<<MinLevelsPerBand<<"  BandWidth = "<<BandWidth<<std::endl;
  out<<" Emin_bands = "<<Emin_bands<<"  Emax_bands = "<<Emax_bands<<std::endl;
  out<<" NLevels = "<<NLevels<<"   NKnownLevels = "<<NKnownLevels<<"   NUnknownLevels = "<<NUnknownLevels<<std::endl;
  out<<" BROpt = "<<BROpt<<"   SampleGammaWidths = "<<SampleGammaWidths<<std::endl;
  out<<" PrimaryGammasIntensityNormFactor = "<<PrimaryGammasIntensityNormFactor<<"   PrimaryGammasEcut = "<<PrimaryGammasEcut<<std::endl;
  out<<" KnownLevelsFlag = "<<KnownLevelsFlag<<std::endl;
  out<<" ElectronConversionFlag = "<<ElectronConversionFlag<<std::endl;
  out<<" ###################################################################################### "<<std::endl;

}

void G4NuDEXStatisticalNucleus::PrintKnownLevels(std::ostream &out){

  out<<" ########################################################################################################## "<<std::endl;
  out<<" KNOWN_LEVEL_SCHEME "<<std::endl;
  out<<" NKnownLevels = "<<NKnownLevels<<std::endl;
  char buffer[1000];

  //for(G4int i=0;i<NKnownLevels;i++){
  for(G4int i=0;i<KnownLevelsVectorSize;i++){
    snprintf(buffer,1000,"%3d %10.4g %5g %2d %10.4g %3d %3d",theKnownLevels[i].id+1,theKnownLevels[i].Energy,theKnownLevels[i].spinx2/2.,2*(G4int)theKnownLevels[i].parity-1,theKnownLevels[i].T12,theKnownLevels[i].NGammas,theKnownLevels[i].Ndecays);
    out<<buffer;
    for(G4int j=0;j<theKnownLevels[i].Ndecays;j++){
      snprintf(buffer,1000,"    %10.4g %7s",theKnownLevels[i].decayFraction[j],theKnownLevels[i].decayMode[j].c_str());
      out<<buffer;
    }
    out<<std::endl;
    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      snprintf(buffer,1000,"                                      %4d %10.4g %10.4g %10.4g %10.4g %10.4g %2d",theKnownLevels[i].FinalLevelID[j]+1,theKnownLevels[i].Eg[j],theKnownLevels[i].Pg[j],theKnownLevels[i].Pe[j],theKnownLevels[i].Icc[j],theKnownLevels[i].cumulPtot[j],theKnownLevels[i].multipolarity[j]);
      out<<buffer<<std::endl;
    }
  }
  out<<" ########################################################################################################## "<<std::endl;

}

void G4NuDEXStatisticalNucleus::PrintKnownLevelsInDEGENformat(std::ostream &out){

  out<<" ########################################################################################################## "<<std::endl;
  out<<" KNOWN_LEVES_DEGEN "<<std::endl;
  out<<" NKnownLevels = "<<NKnownLevels<<std::endl;
  char buffer[1000];

  for(G4int i=0;i<NKnownLevels;i++){
    G4double MaxIntens=-100;
    G4double GammaEnergy;
    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      if(theKnownLevels[i].Pg[j]>MaxIntens){MaxIntens=theKnownLevels[i].Pg[j];}
    }
    for(G4int j=0;j<theKnownLevels[i].NGammas;j++){
      //snprintf(buffer,1000,"%10.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f",theKnownLevels[i].Energy*1000.,theKnownLevels[i].spinx2/2.,2.*(G4int)theKnownLevels[i].parity-1,theKnownLevels[i].Eg[j]*1000.,0.,theKnownLevels[i].Pg[j]/MaxIntens*100.,0.,theKnownLevels[i].Icc[j]);
      GammaEnergy=theKnownLevels[i].Energy-theKnownLevels[theKnownLevels[i].FinalLevelID[j]].Energy;
      snprintf(buffer,1000,"%10.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f",theKnownLevels[i].Energy*1000.,theKnownLevels[i].spinx2/2.,2.*(G4int)theKnownLevels[i].parity-1,GammaEnergy*1000.,0.,theKnownLevels[i].Pg[j]/MaxIntens*100.,0.,theKnownLevels[i].Icc[j]);
      out<<buffer<<std::endl;
    }
  }  
  out<<" ########################################################################################################## "<<std::endl;

}

void G4NuDEXStatisticalNucleus::PrintLevelDensity(std::ostream &out){

  if(theLD==0){return;}

  G4double Emin=0;
  G4double Emax=E_unk_max;
  G4int np=100;

  out<<" ###################################################################################### "<<std::endl;
  out<<" LEVELDENSITY"<<std::endl;
  G4double ene,exp=0;
  G4double *ld=new G4double[maxspinx2+2];
  G4bool *WriteThisSpin=new G4bool[maxspinx2+1];

  for(G4int spx2=0;spx2<=maxspinx2;spx2++){
    WriteThisSpin[spx2]=true;
    if(((A_Int+spx2)%2)!=0){
      WriteThisSpin[spx2]=false;
    }
  }

  out<<np<<"  "<<Emin<<"  "<<Emax<<"  "<<Ecrit<<"  "<<maxspinx2<<std::endl;
  out<<"ENE   EXP   TOT   SUM(J)";
  for(G4int spx2=0;spx2<=maxspinx2;spx2++){
    if(WriteThisSpin[spx2]){out<<"   J="<<spx2/2.;}
  }
  out<<std::endl;

  for(G4int i=0;i<np;i++){
    ene=Emin+(Emax-Emin)*i/(G4double)(np-1);
    exp=0;
    for(G4int j=0;j<NLevels;j++){if(theLevels[j].Energy<ene){exp+=theLevels[j].NLevels;}}
    out<<ene<<"  "<<exp<<"  ";
    ld[maxspinx2+1]=0;
    for(G4int spx2=0;spx2<=maxspinx2;spx2++){
      ld[spx2]=2*theLD->GetLevelDensity(ene,spx2/2.,true);
      ld[maxspinx2+1]+=ld[spx2];
    }
    //out<<ld[maxspinx2+1];
    out<<theLD->GetLevelDensity(ene,0,true,true)<<"  "<<ld[maxspinx2+1];
    for(G4int spx2=0;spx2<=maxspinx2;spx2++){
      if(WriteThisSpin[spx2]){out<<"   "<<ld[spx2];}
    }
    out<<std::endl;
  }
  out<<" ###################################################################################### "<<std::endl;

  delete [] ld;
  delete [] WriteThisSpin;
}

void G4NuDEXStatisticalNucleus::PrintLevelSchemeInDEGENformat(const char* fname,G4int MaxLevelID){

  std::ofstream out(fname);
  char buffer[1000];
  for(G4int i=0;i<NLevels;i++){
    if(theLevels[i].Energy>Ecrit && (MaxLevelID>0 && i<=MaxLevelID)){
      snprintf(buffer,1000,"%13.5f %17.8f %17.8f ",theLevels[i].Energy*1000.,theLevels[i].spinx2/2.,2.*(G4int)theLevels[i].parity-1);
      out<<buffer<<std::endl;
    }
  }
  out.close();

}
  
void G4NuDEXStatisticalNucleus::PrintLevelScheme(std::ostream &out){
  out<<" ###################################################################################### "<<std::endl;
  out<<" LEVELSCHEME"<<std::endl;
  for(G4int i=0;i<NLevels;i++){
    out<<i<<"  "<<theLevels[i].Energy<<"  "<<theLevels[i].spinx2/2.<<"  "<<theLevels[i].parity<<"  "<<theLevels[i].KnownLevelID<<"  "<<theLevels[i].NLevels<<"  "<<theLevels[i].Width<<"  "<<theLevels[i].seed<<std::endl;
  }
  out<<" ###################################################################################### "<<std::endl;
}

void G4NuDEXStatisticalNucleus::PrintThermalPrimaryTransitions(std::ostream &out){

  out<<" #################################################### "<<std::endl;
  out<<" THERMAL PRIMARY TRANSITIONS"<<std::endl;
  out<<" "<<NLevelsBelowThermalCaptureLevel<<std::endl;
  if(theThermalCaptureLevelCumulBR!=0){
    out<<" "<<0<<"  "<<theLevels[0].Energy<<"  "<<Sn-theLevels[0].Energy<<"  "<<theThermalCaptureLevelCumulBR[0]<<std::endl;
    for(G4int i=1;i<NLevelsBelowThermalCaptureLevel;i++){
      out<<" "<<i<<"  "<<theLevels[i].Energy<<"  "<<Sn-theLevels[i].Energy<<"  "<<theThermalCaptureLevelCumulBR[i]-theThermalCaptureLevelCumulBR[i-1]<<std::endl;
    }
  }
  out<<" #################################################### "<<std::endl;

  G4double ThresholdIntensity=0.01;
  out<<" #################################################### "<<std::endl;
  out<<" STRONGEST THERMAL PRIMARY TRANSITIONS"<<std::endl;
  out<<" "<<NLevelsBelowThermalCaptureLevel<<std::endl;
  if(theThermalCaptureLevelCumulBR!=0){
    if(theThermalCaptureLevelCumulBR[0]>ThresholdIntensity){out<<" "<<0<<"  "<<theLevels[0].Energy<<"  "<<Sn-theLevels[0].Energy<<"  "<<theThermalCaptureLevelCumulBR[0]<<std::endl;}
    for(G4int i=1;i<NLevelsBelowThermalCaptureLevel;i++){
      if(theThermalCaptureLevelCumulBR[i]-theThermalCaptureLevelCumulBR[i-1]>ThresholdIntensity){out<<" "<<i<<"  "<<theLevels[i].Energy<<"  "<<Sn-theLevels[i].Energy<<"  "<<theThermalCaptureLevelCumulBR[i]-theThermalCaptureLevelCumulBR[i-1]<<std::endl;}
    }
  }
  out<<" #################################################### "<<std::endl;
}

void G4NuDEXStatisticalNucleus::PrintTotalCumulBR(G4int i_level,std::ostream &out){

  if(TotalCumulBR[i_level]!=0){
    out<<" #################################################### "<<std::endl;
    out<<" CUMULBR FROM LEVEL "<<i_level<<" with ENERGY "<<theLevels[i_level].Energy<<std::endl;
    for(G4int i=0;i<i_level;i++){
      out<<theLevels[i].Energy<<"  "<<theLevels[i].spinx2/2.<<"  "<<theLevels[i].parity<<"  "<<TotalCumulBR[i_level][i]<<std::endl;
    }
    out<<" #################################################### "<<std::endl;
  }

}

void G4NuDEXStatisticalNucleus::PrintBR(G4int i_level,G4double MaxExcEneToPrint_MeV,std::ostream &out){

  if(TotalCumulBR[i_level]!=0){
    out<<" #################################################### "<<std::endl;
    out<<" BR FROM LEVEL "<<i_level<<" with ENERGY "<<theLevels[i_level].Energy<<std::endl;
    for(G4int i=0;i<i_level;i++){
      if(theLevels[i].Energy<MaxExcEneToPrint_MeV || MaxExcEneToPrint_MeV<0){
	if(i==0){
	  out<<theLevels[i].Energy<<"  "<<theLevels[i].spinx2/2.<<"  "<<theLevels[i].parity<<"  "<<TotalCumulBR[i_level][i]<<std::endl;
	}
	else{
	  out<<theLevels[i].Energy<<"  "<<theLevels[i].spinx2/2.<<"  "<<theLevels[i].parity<<"  "<<TotalCumulBR[i_level][i]-TotalCumulBR[i_level][i-1]<<std::endl;
	}
      }
    }
    out<<" #################################################### "<<std::endl;
  }
  
  
}


void G4NuDEXStatisticalNucleus::PrintPSF(std::ostream &out){

  thePSF->PrintPSFParameters(out);

  G4int NVals=400;
  G4int nEnePSF=(G4int)Sn+1; //number of excitation energies where the PSF are evaluated
  G4double EnePSF[200];
  G4double Emin=0;
  G4double Emax=10;
  G4double xval,e1,m1,e2;

  out<<" #################################################### "<<std::endl;
  out<<" PSF"<<std::endl;
  out<<" "<<NVals<<"  "<<Emin<<"  "<<Emax<<"  "<<nEnePSF<<std::endl;
  EnePSF[0]=Sn;
  for(G4int i=1;i<nEnePSF;i++){
    EnePSF[i]=i;
  }
  for(G4int i=0;i<nEnePSF;i++){
    out<<"  "<<EnePSF[i];
  }
  out<<std::endl;
  char word[1000];
  out<<"    E          E1        M1        E2 "<<std::endl;
  for(G4int i=0;i<nEnePSF;i++){
    for(G4int j=0;j<NVals;j++){
      xval=Emin+(Emax-Emin)*j/(NVals-1.);
      if(xval==0){xval=1.e-6;}
      e1=thePSF->GetE1(xval,EnePSF[i]);
      m1=thePSF->GetM1(xval,EnePSF[i]);
      e2=thePSF->GetE2(xval,EnePSF[i]);
      snprintf(word,1000," %10.4E %10.4E %10.4E %10.4E",xval,e1,m1,e2);
      out<<word<<std::endl;
    }
  }
  out<<" #################################################### "<<std::endl;

}


void G4NuDEXStatisticalNucleus::PrintICC(std::ostream &out){

  theICC->PrintICC(out);

}

void G4NuDEXStatisticalNucleus::PrintAll(std::ostream &out){

  PrintParameters(out);
  PrintKnownLevels(out);
  PrintLevelDensity(out);
  PrintLevelScheme(out);
  PrintThermalPrimaryTransitions(out);
  PrintPSF(out);
  PrintICC(out);

}
  


void G4NuDEXStatisticalNucleus::PrintInput01(std::ostream &out){

  out<<"LEVELDENSITYTYPE "<<LevelDensityType<<std::endl;
  out<<"MAXSPIN "<<maxspinx2/2.<<std::endl;
  out<<"MINLEVELSPERBAND "<<MinLevelsPerBand<<std::endl;
  out<<"BANDWIDTH_MEV "<<BandWidth<<std::endl;
  out<<"MAXEXCENERGY_MEV "<<MaxExcEnergy<<std::endl;
  out<<"ECRIT_MEV "<<Ecrit<<std::endl;
  out<<"KNOWNLEVELSFLAG "<<KnownLevelsFlag<<std::endl;
  out<<std::endl;
  out<<"PSF_FLAG "<<PSFflag<<std::endl;
  out<<"BROPTION "<<BROpt<<std::endl;
  out<<"SAMPLEGAMMAWIDTHS "<<SampleGammaWidths<<std::endl;
  out<<std::endl;
  out<<"SEED1 "<<seed1<<std::endl;
  out<<"SEED2 "<<seed2<<std::endl;
  out<<"SEED3 "<<seed3<<std::endl;
  out<<std::endl;
  out<<"ELECTRONCONVERSIONFLAG "<<ElectronConversionFlag<<std::endl;
  out<<"PRIMARYTHCAPGAMNORM "<<PrimaryGammasIntensityNormFactor<<std::endl;
  out<<"PRIMARYGAMMASECUT "<<PrimaryGammasEcut<<std::endl;
  out<<std::endl;
  theLD->PrintParametersInInputFileFormat(out);
  thePSF->PrintPSFParametersInInputFileFormat(out);
  out<<std::endl;
  out<<"END"<<std::endl;

}

G4int ComparisonLevels(const void* va, const void* vb)
{
 Level* a, *b;
 a = (Level*) va;
 b = (Level*) vb;
 if( a->Energy == b->Energy ) return 0;
 return( ( a->Energy ) > ( b->Energy ) ) ? 1:-1;
}

void CopyLevel(Level* a,Level* b){
  b->Energy=a->Energy;
  b->spinx2=a->spinx2;
  b->parity=a->parity;
  b->seed=a->seed;
  b->KnownLevelID=a->KnownLevelID;
  b->NLevels=a->NLevels;
  b->Width=a->Width;
}


void CopyLevel(KnownLevel* a,Level* b){
  b->Energy=a->Energy;
  b->spinx2=a->spinx2;
  b->parity=a->parity;
  b->seed=0;
  b->KnownLevelID=-1;
  b->NLevels=1;
  b->Width=0;
}




