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



#include "G4NuDEXRandom.hh"
#include "G4NuDEXLevelDensity.hh"
#include "G4NuDEXPSF.hh"


G4NuDEXPSF::G4NuDEXPSF(G4int aZ,G4int aA){
  Z_Int=aZ;
  A_Int=aA;
  nR_E1=0; nR_M1=0; nR_E2=0;
  x_E1=0; y_E1=0;
  x_M1=0; y_M1=0;
  x_E2=0; y_E2=0;
  E1_normFac=-1; M1_normFac=-1; E2_normFac=-1;
  NormEmin=0; NormEmax=6; //Integral between 0 and 6 MeV
  ScaleFactor_E1=1;
  ScaleFactor_M1=1;
  ScaleFactor_E2=1;
}

G4NuDEXPSF::~G4NuDEXPSF(){
  if(x_E1!=0){delete [] x_E1;}
  if(y_E1!=0){delete [] y_E1;}
  if(x_M1!=0){delete [] x_M1;}
  if(y_M1!=0){delete [] y_M1;}
  if(x_E2!=0){delete [] x_E2;}
  if(y_E2!=0){delete [] y_E2;}
}


//If inputfname!=0 then we take the PSF data from the inputfname instead of the dirname
G4int G4NuDEXPSF::Init(const char* dirname,G4NuDEXLevelDensity* aLD,const char* inputfname,const char* defaultinputfname,G4int PSFflag){

  theLD=aLD;

  //Three options: very detailed model, if not --> gdr-parameters&errors-exp-MLO.dat (RIPL-3), if not --> theorethical values

  char fname[500];
  G4bool IsDone=false;

  //input:
  if(inputfname!=0){
    IsDone=TakePSFFromInputFile(inputfname);
    if(IsDone){return 0;}
  }

  //default input:
  if(defaultinputfname!=0){
    IsDone=TakePSFFromInputFile(defaultinputfname);
    if(IsDone){return 0;}
  }

  //Detailed model
  snprintf(fname,500,"%s/PSF/PSF_param.dat",dirname);
  IsDone=TakePSFFromDetailedParFile(fname);
  if(IsDone){return 0;}

  //IAEA - 2019 values:
  if(PSFflag==0){
    snprintf(fname,500,"%s/PSF/CRP_IAEA_SMLO_E1_v01.dat",dirname);
    IsDone=TakePSFFromIAEA01(fname);
    if(IsDone){return 0;}
  }
  else if(PSFflag!=1){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  
  //RIPL-MLO values:
  snprintf(fname,500,"%s/PSF/gdr-parameters&errors-exp-MLO.dat",dirname);
  IsDone=TakePSFFromRIPL01(fname);
  if(IsDone){return 0;}

  //RIPL-Theorethical values:
  snprintf(fname,500,"%s/PSF/gdr-parameters-theor.dat",dirname);
  IsDone=TakePSFFromRIPL02(fname);
  if(IsDone){return 0;}

  //Theorethical values:
  // E1 for spherical nucleus:
  nR_E1=0;
  PSFType_E1[nR_E1]=2;
  //G4double a=31.2,b=20.6,c=0.026,d=1.05; //SLO-old (RIPL-2)
  //G4double a=27.47,b=22.063,c=0.0277,d=1.222;//SLO (RIPL-3)
  G4double a=28.69,b=21.731,c=0.0285,d=1.267;//MLO (RIPL-3)

  E_E1[nR_E1]=a*std::pow(A_Int,-1./3.)+b*std::pow(A_Int,-1./6.);
  G_E1[nR_E1]=c*std::pow(E_E1[nR_E1],1.9);
  s_E1[nR_E1]=120/3.141592*d*(A_Int-Z_Int)*Z_Int/(G4double)A_Int/G_E1[nR_E1];
  nR_E1++;
  GenerateM1AndE2FromE1();

  return 0;
}

void  G4NuDEXPSF::GenerateM1AndE2FromE1(){

  //M1:
  nR_M1=0;
  E_M1[nR_M1]=41*std::pow(A_Int,-1./3.);
  G_M1[nR_M1]=4;
  s_M1[nR_M1]=1;
  PSFType_M1[nR_M1]=0;
  nR_M1++;

  //f(E1)/f(M1) = 0.0588*A**0.878    at +-7 MeV
  G4double fE1=GetE1(7,7);
  G4double fM1=GetM1(7,7);
  s_M1[0]=fE1/0.0588/std::pow(A_Int,0.878)/fM1;


  //E2:
  nR_E2=0;
  E_E2[nR_E2]=63*std::pow(A_Int,-1./3.);
  G_E2[nR_E2]=6.11-0.021*A_Int;
  s_E2[nR_E2]=0.00014*Z_Int*Z_Int*E_E2[nR_E2]/std::pow(A_Int,1./3)/G_E2[nR_E2];
  PSFType_E2[nR_E2]=0;
  nR_E2++;

}

G4bool G4NuDEXPSF::TakePSFFromRIPL02(const char* fname){

  G4bool result=false;
  G4int aA,aZ;
  std::ifstream in(fname);
  char dum[200];

  for(G4int i=0;i<4;i++){in.ignore(10000,'\n');}
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      result=true;
      in>>dum>>dum;
      nR_E1=2;
      in>>E_E1[0]>>G_E1[0]>>E_E1[1]>>G_E1[1];
      PSFType_E1[0]=2; PSFType_E1[1]=2; //SMLO 

      G4double a=28.69,b=21.731,c=0.0285,d=1.267;//MLO
      G4double E_E1_0=a*std::pow(A_Int,-1./3.)+b*std::pow(A_Int,-1./6.);
      G4double G_E1_0=c*std::pow(E_E1_0,1.9);
      G4double s_E1_0=120/3.141592*d*(A_Int-Z_Int)*Z_Int/(G4double)A_Int/G_E1_0;
      s_E1[0]=s_E1_0/3.;
      s_E1[1]=2.*s_E1_0/3.;

      break;
    }
    in.ignore(10000,'\n');
  }
  in.close();
  if(result){GenerateM1AndE2FromE1();}

  return result;
}


G4bool G4NuDEXPSF::TakePSFFromRIPL01(const char* fname){

  G4bool result=false;
  G4int aA,aZ;
  std::ifstream in(fname);
  char dum[200];

  for(G4int i=0;i<7;i++){in.ignore(10000,'\n');}
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      result=true;
      in>>dum>>dum;
      nR_E1=0;
      in>>E_E1[nR_E1]>>s_E1[nR_E1]>>G_E1[nR_E1];
      PSFType_E1[nR_E1]=2; //SMLO
      nR_E1++;
      //sometimes there is a second resonance:
      in>>E_E1[nR_E1]>>dum>>G_E1[nR_E1];
      if(dum[0]!='-'){ //there is a second resonance
	s_E1[nR_E1]=std::atof(dum);
	PSFType_E1[nR_E1]=2;
	nR_E1++;
      }
      break;
    }
    in.ignore(10000,'\n');
  }
  in.close();
  if(result){GenerateM1AndE2FromE1();}

  return result;
}


G4bool G4NuDEXPSF::TakePSFFromIAEA01(const char* fname){

  G4bool result=false;
  G4int aA,aZ;
  char dum[200];
  G4double beta=0;
  std::ifstream in(fname);
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      result=true;
      nR_E1=0;
      in>>dum>>dum>>E_E1[nR_E1]>>dum>>dum>>G_E1[nR_E1]>>dum>>dum>>s_E1[nR_E1];
      PSFType_E1[nR_E1]=11;
      nR_E1++;
      in>>dum;
      if(std::string(dum)==std::string("beta=")){
	in>>beta;
	break;
      }
      else if(std::string(dum)==std::string("Er2")){
	in>>dum>>E_E1[nR_E1]>>dum>>dum>>G_E1[nR_E1]>>dum>>dum>>s_E1[nR_E1]>>dum>>beta;
	PSFType_E1[nR_E1]=11;
	nR_E1++;

      }
      else{
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      break;
    }
    in.ignore(10000,'\n');
  }
  if(!result){
    return result;
  }

  //---------------------------------------------------
  //Now M1 (https://doi.org/10.1140/epja/i2019-12840-1)
  nR_M1=0;
  //Spin-flip:
  PSFType_M1[nR_M1]=0;
  E_M1[nR_M1]=18.0*std::pow(A_Int,-1./6.);
  G_M1[nR_M1]=4;
  s_M1[nR_M1]=0.03*std::pow(A_Int,5./6.);
  nR_M1++;
  //Scissors-mode:
  PSFType_M1[nR_M1]=0;
  E_M1[nR_M1]=5.0*std::pow(A_Int,-1./10.);
  G_M1[nR_M1]=1.5;
  s_M1[nR_M1]=0.02*std::fabs(beta)*std::pow(A_Int,9./10.);
  nR_M1++;
  //upbend:
  PSFType_M1[nR_M1]=21;
  E_M1[nR_M1]=0.4035*std::exp(-6.0*std::fabs(beta));
  G_M1[nR_M1]=0.8;
  s_M1[nR_M1]=0;
  nR_M1++;
  //---------------------------------------------------

  //---------------------------------------------------
  //E2 same as in the old RIPL recommendations:
  nR_E2=0;
  E_E2[nR_E2]=63*std::pow(A_Int,-1./3.);
  G_E2[nR_E2]=6.11-0.021*A_Int;
  s_E2[nR_E2]=0.00014*Z_Int*Z_Int*E_E2[nR_E2]/std::pow(A_Int,1./3)/G_E2[nR_E2];
  PSFType_E2[nR_E2]=0;
  nR_E2++;
  //---------------------------------------------------
  
  return result;
}



G4bool G4NuDEXPSF::TakePSFFromInputFile(const char* fname){

  G4bool result=false;
  char word[1000];
  std::ifstream in(fname);
  while(in>>word){
    if(word[0]=='#'){in.ignore(10000,'\n');}
    if(std::string(word)==std::string("END")){break;}
    if(std::string(word)==std::string("PSF")){
      result=true;
      in>>nR_E1;
      for(G4int i=0;i<nR_E1;i++){
	in>>PSFType_E1[i]>>E_E1[i]>>G_E1[i]>>s_E1[i];
	if(PSFType_E1[i]==7){in>>p1_E1[i];}
	if(PSFType_E1[i]==8){in>>p1_E1[i]>>p2_E1[i];}
	if(PSFType_E1[i]==9){in>>p1_E1[i]>>p2_E1[i];}
	if(PSFType_E1[i]==10){in>>p1_E1[i]>>p2_E1[i]>>p3_E1[i];}
	if(PSFType_E1[i]==40 || PSFType_E1[i]==41){ //only one pointwise function is allowed
	  if(x_E1!=0){NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");}
	  in>>np_E1;
	  x_E1=new G4double[np_E1]; y_E1=new G4double[np_E1]; 
	  for(G4int j=0;j<np_E1;j++){in>>x_E1[j]>>y_E1[j];}
	  in>>E1_normFac;
	}
      }
      in>>nR_M1;
      for(G4int i=0;i<nR_M1;i++){
	in>>PSFType_M1[i]>>E_M1[i]>>G_M1[i]>>s_M1[i];
	if(PSFType_M1[i]==7){in>>p1_M1[i];}
	if(PSFType_M1[i]==8){in>>p1_M1[i]>>p2_M1[i];}
	if(PSFType_M1[i]==9){in>>p1_M1[i]>>p2_M1[i];}
	if(PSFType_M1[i]==10){in>>p1_M1[i]>>p2_M1[i]>>p3_M1[i];}
	if(PSFType_M1[i]==40 || PSFType_M1[i]==41){//only one pointwise function is allowed
	  if(x_M1!=0){NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");}
	  in>>np_M1;
	  x_M1=new G4double[np_M1]; y_M1=new G4double[np_M1]; 
	  for(G4int j=0;j<np_M1;j++){in>>x_M1[j]>>y_M1[j];}
	  in>>M1_normFac;
	}
      }
      in>>nR_E2;
      for(G4int i=0;i<nR_E2;i++){
	in>>PSFType_E2[i]>>E_E2[i]>>G_E2[i]>>s_E2[i];
	if(PSFType_E2[i]==7){in>>p1_E2[i];}
	if(PSFType_E2[i]==8){in>>p1_E2[i]>>p2_E2[i];}
	if(PSFType_E2[i]==9){in>>p1_E2[i]>>p2_E2[i];}
	if(PSFType_E2[i]==10){in>>p1_E2[i]>>p2_E2[i]>>p3_E2[i];}
	if(PSFType_E2[i]==40 || PSFType_E2[i]==41){//only one pointwise function is allowed
	  if(x_E2!=0){NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");}
	  in>>np_E2;
	  x_E2=new G4double[np_E2]; y_E2=new G4double[np_E2]; 
	  for(G4int j=0;j<np_E2;j++){in>>x_E2[j]>>y_E2[j];}
	  in>>E2_normFac;
	}
      }
      break;
    }
  }
  
  Renormalize(); // if XX_normFac>0 --> renormalization of the PSF

  return result;
}


void G4NuDEXPSF::Renormalize(){

  G4int npIntegral=1000;
  G4double Integral=0,x_eval,y_eval;
  G4double binWidth=(NormEmax-NormEmin)/npIntegral;

  //-------------------------------------------------
  if(E1_normFac>0){
    Integral=0;
    for(G4int i=0;i<npIntegral;i++){
      x_eval=NormEmin+binWidth*(i+0.5);
      y_eval=GetE1(x_eval,NormEmax);
      Integral+=y_eval;
    }
    Integral*=binWidth;
    ScaleFactor_E1=E1_normFac/Integral;
  }
  //-------------------------------------------------
  //-------------------------------------------------
  if(M1_normFac>0){
    Integral=0;
    for(G4int i=0;i<npIntegral;i++){
      x_eval=NormEmin+binWidth*(i+0.5);
      y_eval=GetM1(x_eval,NormEmax);
      Integral+=y_eval;
    }
    Integral*=binWidth;
    ScaleFactor_M1=M1_normFac/Integral;
    //std::cout<<M1_normFac<<"  "<<Integral<<"  "<<ScaleFactor_M1<<std::endl; getchar();
  }
  //-------------------------------------------------
  //-------------------------------------------------
  if(E2_normFac>0){
    Integral=0;
    for(G4int i=0;i<npIntegral;i++){
      x_eval=NormEmin+binWidth*(i+0.5);
      y_eval=GetE2(x_eval,NormEmax);
      Integral+=y_eval;
    }
    Integral*=binWidth;
    ScaleFactor_E2=E2_normFac/Integral;
  }
  //-------------------------------------------------

}



G4bool G4NuDEXPSF::TakePSFFromDetailedParFile(const char* fname){

  G4bool result=false;
  G4int aA,aZ;
  std::ifstream in(fname);
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      result=true;
      in>>nR_E1;
      for(G4int i=0;i<nR_E1;i++){
	in>>PSFType_E1[i]>>E_E1[i]>>G_E1[i]>>s_E1[i];
	if(PSFType_E1[i]==7){in>>p1_E1[i];}
	if(PSFType_E1[i]==8){in>>p1_E1[i]>>p2_E1[i];}
	if(PSFType_E1[i]==9){in>>p1_E1[i]>>p2_E1[i];}
	if(PSFType_E1[i]==10){in>>p1_E1[i]>>p2_E1[i]>>p3_E1[i];}
      }
      in>>nR_M1;
      for(G4int i=0;i<nR_M1;i++){
	in>>PSFType_M1[i]>>E_M1[i]>>G_M1[i]>>s_M1[i];
	if(PSFType_M1[i]==7){in>>p1_M1[i];}
	if(PSFType_M1[i]==8){in>>p1_M1[i]>>p2_M1[i];}
 	if(PSFType_M1[i]==9){in>>p1_M1[i]>>p2_M1[i];}
	if(PSFType_M1[i]==10){in>>p1_M1[i]>>p2_M1[i]>>p3_M1[i];}
     }
      in>>nR_E2;
      for(G4int i=0;i<nR_E2;i++){
	in>>PSFType_E2[i]>>E_E2[i]>>G_E2[i]>>s_E2[i];
 	if(PSFType_E2[i]==7){in>>p1_E2[i];}
	if(PSFType_E2[i]==8){in>>p1_E2[i]>>p2_E2[i];}
	if(PSFType_E2[i]==9){in>>p1_E2[i]>>p2_E2[i];}
	if(PSFType_E2[i]==10){in>>p1_E2[i]>>p2_E2[i]>>p3_E2[i];}
     }
      break;
    }
    in.ignore(10000,'\n');
  }
  in.close();

  return result;
}





G4double G4NuDEXPSF::GetE1(G4double Eg,G4double ExcitationEnergy){

  G4double result=0;
  for(G4int i=0;i<nR_E1;i++){
    if(PSFType_E1[i]==0){
      result+=8.674E-8*SLO(Eg,E_E1[i],G_E1[i],s_E1[i]);
    }
    else if(PSFType_E1[i]==1){
      result+=8.674E-8*EGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==2){
      result+=8.674E-8*SMLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==3){
      result+=8.674E-8*GLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==4){
      result+=8.674E-8*MGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==5){
      result+=8.674E-8*KMF(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==6){
      result+=8.674E-8*GH(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==7){
      result+=8.674E-8*MEGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy,p1_E1[i],p1_E1[i]);
    }
    else if(PSFType_E1[i]==8){
      result+=8.674E-8*MEGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy,p1_E1[i],p2_E1[i]);
    }
    else if(PSFType_E1[i]==9){
      result+=8.674E-8*MEGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy,p1_E1[i],p1_E1[i],p2_E1[i]);
    }
    else if(PSFType_E1[i]==10){
      result+=8.674E-8*MEGLO(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy,p1_E1[i],p2_E1[i],p3_E1[i]);
    }
    else if(PSFType_E1[i]==11){
      result+=8.674E-8*SMLO_v2(Eg,E_E1[i],G_E1[i],s_E1[i],ExcitationEnergy);
    }
    else if(PSFType_E1[i]==20){
      result+=8.674E-8*Gauss(Eg,E_E1[i],G_E1[i],s_E1[i]);
    }
    else if(PSFType_E1[i]==21){
      result+=8.674E-8*Expo(Eg,E_E1[i],G_E1[i]);
    }
    else if(PSFType_E1[i]==40){
      result+=EvaluateFunction(Eg,np_E1,x_E1,y_E1);
    }
    else if(PSFType_E1[i]==41){
      result+=std::pow(10.,EvaluateFunction(Eg,np_E1,x_E1,y_E1));
    }
    else{
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
  }

  if(result!=result){ // nan
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  return result*ScaleFactor_E1;
}

G4double G4NuDEXPSF::GetM1(G4double Eg,G4double ExcitationEnergy){

  G4double result=0;
  for(G4int i=0;i<nR_M1;i++){
    if(PSFType_M1[i]==0){
      result+=8.674E-8*SLO(Eg,E_M1[i],G_M1[i],s_M1[i]);
    }
    else if(PSFType_M1[i]==1){
      result+=8.674E-8*EGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==2){
      result+=8.674E-8*SMLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==3){
      result+=8.674E-8*GLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==4){
      result+=8.674E-8*MGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==5){
      result+=8.674E-8*KMF(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==6){
      result+=8.674E-8*GH(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==7){
      result+=8.674E-8*MEGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy,p1_M1[i],p1_M1[i]);
    }
    else if(PSFType_M1[i]==8){
      result+=8.674E-8*MEGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy,p1_M1[i],p2_M1[i]);
    }
    else if(PSFType_M1[i]==9){
      result+=8.674E-8*MEGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy,p1_M1[i],p1_M1[i],p2_M1[i]);
    }
    else if(PSFType_M1[i]==10){
      result+=8.674E-8*MEGLO(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy,p1_M1[i],p2_M1[i],p3_M1[i]);
    }
    else if(PSFType_M1[i]==11){
      result+=8.674E-8*SMLO_v2(Eg,E_M1[i],G_M1[i],s_M1[i],ExcitationEnergy);
    }
    else if(PSFType_M1[i]==20){
      result+=8.674E-8*Gauss(Eg,E_M1[i],G_M1[i],s_M1[i]);
    }
    else if(PSFType_M1[i]==21){
      result+=8.674E-8*Expo(Eg,E_M1[i],G_M1[i]);
    }
    else if(PSFType_M1[i]==40){
      result+=EvaluateFunction(Eg,np_M1,x_M1,y_M1);
    }
    else if(PSFType_M1[i]==41){
      result+=std::pow(10.,EvaluateFunction(Eg,np_M1,x_M1,y_M1));
    }
    else{
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
  }

  if(result!=result){ // nan
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  return result*ScaleFactor_M1;
}

G4double G4NuDEXPSF::GetE2(G4double Eg,G4double ExcitationEnergy){

  G4double result=0;
  for(G4int i=0;i<nR_E2;i++){
    if(PSFType_E2[i]==0){
      result+=5.22E-8*SLO(Eg,E_E2[i],G_E2[i],s_E2[i]);
    }
    else if(PSFType_E2[i]==1){
      result+=5.22E-8*EGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==2){
      result+=5.22E-8*SMLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==3){
      result+=5.22E-8*GLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==4){
      result+=5.22E-8*MGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==5){
      result+=5.22E-8*KMF(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==6){
      result+=5.22E-8*GH(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==7){
      result+=5.22E-8*MEGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy,p1_E2[i],p1_E2[i]);
    }
    else if(PSFType_E2[i]==8){
      result+=5.22E-8*MEGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy,p1_E2[i],p2_E2[i]);
    }
    else if(PSFType_E2[i]==9){
      result+=5.22E-8*MEGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy,p1_E2[i],p1_E2[i],p2_E2[i]);
    }
    else if(PSFType_E2[i]==10){
      result+=5.22E-8*MEGLO(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy,p1_E2[i],p2_E2[i],p3_E2[i]);
    }
    else if(PSFType_E2[i]==11){
      result+=5.22E-8*SMLO_v2(Eg,E_E2[i],G_E2[i],s_E2[i],ExcitationEnergy);
    }
    else if(PSFType_E2[i]==20){
      result+=5.22E-8*Gauss(Eg,E_E2[i],G_E2[i],s_E2[i]);
    }
    else if(PSFType_E2[i]==21){
      result+=5.22E-8*Expo(Eg,E_E2[i],G_E2[i]);
    }
    else if(PSFType_E2[i]==40){
      result+=EvaluateFunction(Eg,np_E2,x_E2,y_E2);
    }
    else if(PSFType_E2[i]==41){
      result+=std::pow(10.,EvaluateFunction(Eg,np_E2,x_E2,y_E2));
    }
    else{
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
  }

  if(result!=result){ // nan
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  return result*ScaleFactor_E2;
}


//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************

//Defined as in RIPL-3 when possible. Some of them come from other references:
//http://dx.doi.org/10.1103/PhysRevC.88.034317

G4double G4NuDEXPSF::SLO(G4double Eg,G4double Er,G4double Gr,G4double sr){

  return sr*Gr*Eg*Gr/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gr*Gr);

}

//Kadmenskij-Markushev-Furman model (KMF) --> not well described in RIPL-3 manual, taken from another document
G4double G4NuDEXPSF::KMF(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double Tf=0;
  if(theLD!=0){
    Tf=theLD->GetNucleusTemperature(ExcitationEnergy-Eg);
  }
  G4double Gc=Gr/Er/Er*(Eg*Eg+4*3.141592*3.141592*Tf*Tf);

  if(Eg==Er){return 0;}
  
  return 0.7*Er*Gr*sr*Gc/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er));
}


G4double G4NuDEXPSF::EGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double result=EGLO_GLO_MGLO(Eg,Er,Gr,sr,ExcitationEnergy,0);

  return result;  
}

G4double G4NuDEXPSF::GLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double result=EGLO_GLO_MGLO(Eg,Er,Gr,sr,ExcitationEnergy,1);

  return result;  
}

G4double G4NuDEXPSF::MGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double result=EGLO_GLO_MGLO(Eg,Er,Gr,sr,ExcitationEnergy,2);

  return result;  
}

//Hybrid model
G4double G4NuDEXPSF::GH(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double Tf=0;
  if(theLD!=0){
    Tf=theLD->GetNucleusTemperature(ExcitationEnergy-Eg);
  }

  G4double Gamma_h=0.63*Gr/Eg/Er*(Eg*Eg+4*3.141592*3.141592*Tf*Tf);
    
  return sr*Gr*Eg*Gamma_h/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gr*Gamma_h);

}

G4double G4NuDEXPSF::SMLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double Tf=0;
  if(theLD!=0){
    Tf=theLD->GetNucleusTemperature(ExcitationEnergy-Eg);
  }
  
  G4double Lambda=1/(1.-std::exp(-Eg/Tf)); 
  G4double Gk_Eg=Gr/Er*ExcitationEnergy;

  return Lambda*sr*Gr*Eg*Gk_Eg/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gk_Eg*Gk_Eg);
}


G4double G4NuDEXPSF::SMLO_v2(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy){

  G4double Tf=0;
  if(Eg<ExcitationEnergy){
    Tf=std::sqrt((ExcitationEnergy-Eg)/(A_Int/10.));
  }
  
  G4double Lambda=1/(1.-std::exp(-Eg/Tf));
  G4double sig_trk=60.*(A_Int-Z_Int)*Z_Int/(G4double)A_Int;
  G4double Gk_Eg=Gr/Er*(Eg+4*3.141592*3.141592*Tf*Tf/Er);

  return Lambda*sig_trk*2./3.141592*sr*Eg*Gk_Eg/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gk_Eg*Gk_Eg);
}


G4double G4NuDEXPSF::Gauss(G4double Eg,G4double Er,G4double sigma,G4double Area){

  return Area*(1./(sigma*std::sqrt(2.*3.141592)))*std::exp(-0.5*std::pow((Eg-Er)/sigma,2.));

}

G4double G4NuDEXPSF::Expo(G4double Eg,G4double C,G4double eta){

  return C*std::exp(-eta*Eg);

}


G4double G4NuDEXPSF::MEGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy,G4double k_param1,G4double k_param2,G4double Temp){

  G4double /*Ti=0,*/Tf=0;
  if(Temp>=0){
    //Ti=Temp;
    Tf=Temp;
  }
  else if(theLD!=0){
    //Ti=theLD->GetNucleusTemperature(ExcitationEnergy);
    Tf=theLD->GetNucleusTemperature(ExcitationEnergy-Eg);
  }

  G4double Gk_Eg=Gamma_k(Eg,Er,Gr,Tf,k_param1);
  //G4double Gk_0=Gamma_k(0,Er,Gr,Ti,k_param2);
  G4double Gk_0=Gamma_k(0,Er,Gr,Tf,k_param2); // in most of the references they use just one temperature

  return sr*Gr*(Eg*Gk_Eg/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gk_Eg*Gk_Eg)+0.7*Gk_0/Er/Er/Er);

}




//Ti, Tf  --> initial/final temperature of the nucleus
//Opt = 0,1,2 --> EGLO, GLO, MGLO
G4double G4NuDEXPSF::EGLO_GLO_MGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy,G4int Opt){

  G4double Ti=0,Tf=0;
  if(theLD!=0){
    Ti=theLD->GetNucleusTemperature(ExcitationEnergy);
    Tf=theLD->GetNucleusTemperature(ExcitationEnergy-Eg);
  }
  
  //k_param could be modified according to experimental data.
  //The following expression is just a general recomendation
  //If k_param==1 --> GLO
  G4double k_param=1; 
  if(A_Int>=148){
    k_param=1+0.09*(A_Int-148)*(A_Int-148)*std::exp(-0.18*(A_Int-148));
  }
  G4double result=0;
  if(Opt==0){//EGLO
    result=FlexibleGLOType(Eg,Er,Gr,sr,Tf,k_param,Ti,k_param);
  }
  else if(Opt==1){//GLO --> same as EGLO, but k_param=1
    result=FlexibleGLOType(Eg,Er,Gr,sr,Tf,1,Ti,1);
  }
  else if(Opt==2){//MGLO --> same as EGLO, but k_param2=1
    result=FlexibleGLOType(Eg,Er,Gr,sr,Tf,k_param,Ti,1);
  }
  else{
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  
  return result;  
}


G4double G4NuDEXPSF::FlexibleGLOType(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double Temp1,G4double k_param1,G4double /*Temp2*/,G4double k_param2){

  G4double Gk_Eg=Gamma_k(Eg,Er,Gr,Temp1,k_param1);
  //G4double Gk_0=Gamma_k(0,Er,Gr,Temp2,k_param2);
  G4double Gk_0=Gamma_k(0,Er,Gr,Temp1,k_param2); // in most of the references they use just one temperature

  return sr*Gr*(Eg*Gk_Eg/((Eg*Eg-Er*Er)*(Eg*Eg-Er*Er)+Eg*Eg*Gk_Eg*Gk_Eg)+0.7*Gk_0/Er/Er/Er);

}


G4double G4NuDEXPSF::Gamma_k(G4double Eg,G4double Er,G4double Gr,G4double Temp,G4double k_param){

  G4double eps0_param=4.5;
  G4double Chi=1;
  if(Er>eps0_param){
    Chi=k_param+(1-k_param)*(Eg-eps0_param)/(Er-eps0_param);
  }
  G4double C_coll=Gr/Er/Er*Chi;
  G4double Gamma_k=C_coll*(Eg*Eg+4*3.141592*3.141592*Temp*Temp);

  return Gamma_k;
}




//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************



void G4NuDEXPSF::PrintPSFParameters(std::ostream &out){

  out<<" ###################################################################################### "<<std::endl;
  out<<" PSF_PARAMS"<<std::endl;
  out<<" E1: nRes = "<<nR_E1<<std::endl;
  for(G4int i=0;i<nR_E1;i++){
    out<<"   "<<PSFType_E1[i]<<"  "<<E_E1[i]<<"  "<<G_E1[i]<<"  "<<s_E1[i]<<std::endl;
    if(PSFType_E1[i]==7){out<<"                       "<<p1_E1[i]<<std::endl;}
    if(PSFType_E1[i]==8){out<<"                       "<<p1_E1[i]<<"  "<<p2_E1[i]<<std::endl;}
    if(PSFType_E1[i]==9){out<<"                       "<<p1_E1[i]<<"  "<<p2_E1[i]<<std::endl;}
    if(PSFType_E1[i]==10){out<<"                       "<<p1_E1[i]<<"  "<<p2_E1[i]<<"  "<<p3_E1[i]<<std::endl;}
    if(PSFType_E1[i]==40 || PSFType_E1[i]==41){out<<np_E1; for(G4int j=0;j<np_E1;j++){out<<"  "<<x_E1[j]<<"  "<<y_E1[j];} out<<std::endl;}
  }
  out<<" M1: nRes = "<<nR_M1<<std::endl;
  for(G4int i=0;i<nR_M1;i++){
    out<<"   "<<PSFType_M1[i]<<"  "<<E_M1[i]<<"  "<<G_M1[i]<<"  "<<s_M1[i]<<std::endl;
    if(PSFType_M1[i]==7){out<<"                       "<<p1_M1[i]<<std::endl;}
    if(PSFType_M1[i]==8){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i]<<std::endl;}
    if(PSFType_M1[i]==9){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i]<<std::endl;}
    if(PSFType_M1[i]==10){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i]<<"  "<<p3_M1[i]<<std::endl;}
    if(PSFType_M1[i]==40 || PSFType_M1[i]==41){out<<np_M1; for(G4int j=0;j<np_M1;j++){out<<"  "<<x_M1[j]<<"  "<<y_M1[j];} out<<std::endl;}
  }
  out<<" E2: nRes = "<<nR_E2<<std::endl;
  for(G4int i=0;i<nR_E2;i++){
    out<<"   "<<PSFType_E2[i]<<"  "<<E_E2[i]<<"  "<<G_E2[i]<<"  "<<s_E2[i]<<std::endl;
    if(PSFType_E2[i]==7){out<<"                       "<<p1_E2[i]<<std::endl;}
    if(PSFType_E2[i]==8){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i]<<std::endl;}
    if(PSFType_E2[i]==9){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i]<<std::endl;}
    if(PSFType_E2[i]==10){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i]<<"  "<<p3_E2[i]<<std::endl;}
    if(PSFType_E2[i]==40 || PSFType_E2[i]==41){out<<np_E2; for(G4int j=0;j<np_E2;j++){out<<"  "<<x_E2[j]<<"  "<<y_E2[j];} out<<std::endl;}
  }
  out<<" ###################################################################################### "<<std::endl;

}


void G4NuDEXPSF::PrintPSFParametersInInputFileFormat(std::ostream &out){

  out<<" PSF"<<std::endl;
  out.precision(15);
  out<<nR_E1<<std::endl;
  for(G4int i=0;i<nR_E1;i++){
    out<<"   "<<PSFType_E1[i]<<"  "<<E_E1[i]<<"  "<<G_E1[i]<<"  "<<s_E1[i];
    if(PSFType_E1[i]==7){out<<"   "<<p1_E1[i];}
    if(PSFType_E1[i]==8){out<<"   "<<p1_E1[i]<<"  "<<p2_E1[i];}
    if(PSFType_E1[i]==9){out<<"   "<<p1_E1[i]<<"  "<<p2_E1[i];}
    if(PSFType_E1[i]==10){out<<"  "<<p1_E1[i]<<"  "<<p2_E1[i]<<"  "<<p3_E1[i];}
    if(PSFType_E1[i]==40 || PSFType_E1[i]==41){out<<np_E1; for(G4int j=0;j<np_E1;j++){out<<"  "<<x_E1[j]<<"  "<<y_E1[j];} }
    out<<std::endl;
  }
  out<<nR_M1<<std::endl;
  for(G4int i=0;i<nR_M1;i++){
    out<<"   "<<PSFType_M1[i]<<"  "<<E_M1[i]<<"  "<<G_M1[i]<<"  "<<s_M1[i];
    if(PSFType_M1[i]==7){out<<"                       "<<p1_M1[i];}
    if(PSFType_M1[i]==8){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i];}
    if(PSFType_M1[i]==9){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i];}
    if(PSFType_M1[i]==10){out<<"                       "<<p1_M1[i]<<"  "<<p2_M1[i]<<"  "<<p3_M1[i];}
    if(PSFType_M1[i]==40 || PSFType_M1[i]==41){out<<np_M1; for(G4int j=0;j<np_M1;j++){out<<"  "<<x_M1[j]<<"  "<<y_M1[j];}}
    out<<std::endl;
  }
  out<<nR_E2<<std::endl;
  for(G4int i=0;i<nR_E2;i++){
    out<<"   "<<PSFType_E2[i]<<"  "<<E_E2[i]<<"  "<<G_E2[i]<<"  "<<s_E2[i];
    if(PSFType_E2[i]==7){out<<"                       "<<p1_E2[i];}
    if(PSFType_E2[i]==8){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i];}
    if(PSFType_E2[i]==9){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i];}
    if(PSFType_E2[i]==10){out<<"                       "<<p1_E2[i]<<"  "<<p2_E2[i]<<"  "<<p3_E2[i];}
    if(PSFType_E2[i]==40 || PSFType_E2[i]==41){out<<np_E2; for(G4int j=0;j<np_E2;j++){out<<"  "<<x_E2[j]<<"  "<<y_E2[j];}}
    out<<std::endl;
  }

}

G4double G4NuDEXPSF::EvaluateFunction(G4double xval,G4int np,G4double* x,G4double* y){

  if(xval<x[0]){return y[0];}
  if(xval>x[np-1]){return y[np-1];}

  G4double m,b;
  G4int i_eval=np-1;
  for(G4int i=1;i<np;i++){
    if(x[i]>=xval){
      i_eval=i;
      break;
    }
  }

  m=(y[i_eval]-y[i_eval-1])/(x[i_eval]-x[i_eval-1]);
  b=y[i_eval]-m*x[i_eval];

  return m*xval+b;
}

