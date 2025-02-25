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



G4NuDEXLevelDensity::G4NuDEXLevelDensity(G4int aZ,G4int aA,G4int ldtype){

  Z_Int=aZ;
  A_Int=aA;
  LDType=ldtype;
  if(LDType<0){LDType=DEFAULTLDTYPE;}
  A_mass=A_Int;
  if(LDType!=1 && LDType!=2 && LDType!=3){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  HasData=false;

  Sn=-1; D0=-1; I0=-1000;
  Ed=0; 
  ainf_ldpar=0; gamma_ldpar=0; dW_ldpar=0; Delta_ldpar=0; T_ldpar=0; E0_ldpar=0; Ex_ldpar=0;
}


G4int G4NuDEXLevelDensity::ReadLDParameters(const char* dirname,const char* inputfname,const char* defaultinputfname){

  char fname[100];
  if(LDType==1 || LDType==3){ // Back-Shifted-Fermi-Gas model
    snprintf(fname,100,"%s/LevelDensities/level-densities-bfmeff.dat",dirname);
  }
  else{ //  Constant Temperature
    snprintf(fname,100,"%s/LevelDensities/level-densities-ctmeff.dat",dirname);
  }
  G4double EL=0,EU=0;
  
  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  G4int aZ,aA;
  char word[200];
  in.ignore(10000,'\n');

  //std::cout<<" LDType="<<LDType<<"  "<<fname<<"  "<<Z_Int<<"  "<<A_Int<<std::endl;

  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      if(LDType==1 || LDType==3){
	in>>word>>I0>>Sn>>D0>>word>>word>>EL>>word>>EU>>dW_ldpar>>gamma_ldpar>>ainf_ldpar>>word>>Delta_ldpar;
	Ex_ldpar=0; E0_ldpar=0; T_ldpar=0;
	Ed=(EL+EU)/2.;
      }
      else if(LDType==2){
	in>>word>>I0>>Sn>>D0>>word>>word>>EL>>word>>EU>>dW_ldpar>>gamma_ldpar>>ainf_ldpar>>word>>Delta_ldpar>>Ex_ldpar>>E0_ldpar>>T_ldpar;
	Ed=(EL+EU)/2.;
      } 
      else{
	NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
      }
      if(in.good()){
	HasData=true;
	break;
      }
    }
    in.ignore(10000,'\n');
  }
  in.close();


  //Re-write some parameters if inputfname!=0
  if(defaultinputfname!=0){
    SearchLDParametersInInputFile(defaultinputfname);
  }
  if(inputfname!=0){
    SearchLDParametersInInputFile(inputfname);
  }

  if(!HasData){ //no data
    G4int check=CalculateLDParameters_BSFG(dirname);
    if(check==0){
      HasData=true;
      if(LDType==2){
	LDType=1;
	std::cout<<" ##### WARNING: level density model for ZA="<<Z_Int*1000+A_Int<<" changed to Back-Shifted-Fermi-Gas model #####"<<std::endl;
      }
    }
  }
  
  if(HasData){return 0;} 

  
  //else, some problem reading ...
  return -1;
}


G4int G4NuDEXLevelDensity::CalculateLDParameters_BSFG(const char* dirname){

  //Eq. 61 of RIPL-3:
  G4double alpha=0.0722396; //MeV^{-1}
  G4double beta= 0.195267; //MeV^{-1}
  G4double gamma0=0.410289; //MeV^{-1}
  G4double delta=0.173015; //MeV

  //Delta_ldpar: Eq. 50 of RIPL-3:
  G4double n_par=0;
  if((Z_Int%2)==1 && ((A_Int-Z_Int)%2)==1){n_par=-1;} //odd-odd (impar-impar)
  if((Z_Int%2)==0 && ((A_Int-Z_Int)%2)==0){n_par=1;} //even-even (par-par)
  Delta_ldpar=n_par*12/std::sqrt(A_mass)+delta;

  //ainf_ldpar: Eq. 52 of RIPL-3:
  ainf_ldpar=alpha*A_Int+beta*std::pow(A_mass,2./3.);

  //gamma_ldpar: Eq. 53 of RIPL-3:
  gamma_ldpar=gamma0/std::pow(A_mass,1./3.);

  //dW_ldpar --> from data file
  char fname[100];
  snprintf(fname,100,"%s/LevelDensities/shellcor-ms.dat",dirname);
  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ######## Error opening file "<<fname<<" ########"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  G4int aZ,aA;
  char word[200];
  in.ignore(10000,'\n');
  in.ignore(10000,'\n');
  in.ignore(10000,'\n');
  in.ignore(10000,'\n');
  while(in>>aZ>>aA){
    if(aZ==Z_Int && aA==A_Int){
      in>>word>>dW_ldpar;
      if(in.good()){break;}
    }
    in.ignore(10000,'\n');
  }
  if(!in.good()){//no data found
    return -1;
  }
  in.close();

  Ex_ldpar=0; E0_ldpar=0; T_ldpar=0;
  Ed=0;


  return 0;
}


G4int G4NuDEXLevelDensity::SearchLDParametersInInputFile(const char* inputfname){

  if(inputfname!=0){
    std::ifstream in2(inputfname);
    if(!in2.good()){
      std::cout<<" ############## Error opening "<<inputfname<<" ##############"<<std::endl;
      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
    }
    std::string word_tmp;
    while(in2>>word_tmp){
      if(word_tmp[0]=='#'){in2.ignore(10000,'\n');}
      if(word_tmp==std::string("END")){break;}
      if(word_tmp==std::string("LDPARAMETERS")){
	in2>>LDType;
	if(LDType==1){
	  in2>>dW_ldpar>>gamma_ldpar>>ainf_ldpar>>Delta_ldpar;
	}
	else if(LDType==2){
	  in2>>dW_ldpar>>gamma_ldpar>>ainf_ldpar>>Delta_ldpar>>Ex_ldpar>>E0_ldpar>>T_ldpar;
	}
	else if(LDType==3){
	  in2>>ainf_ldpar>>Delta_ldpar;
	}
	else{
	  std::cout<<" ############## Error: Unknown LDType="<<LDType<<" in "<<inputfname<<" ##############"<<std::endl;
	  NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	}
	if(!in2.good()){
	  std::cout<<" ############## Error reading "<<inputfname<<" ##############"<<std::endl;
	  NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	}
	HasData=true;
	break;
      }
    }
    in2.close();
  }
  
  return 0;
}

void G4NuDEXLevelDensity::PrintParametersInInputFileFormat(std::ostream &out){

  out<<"LDPARAMETERS"<<std::endl;
  out<<LDType<<std::endl;
  G4long oldprc = out.precision(15);
  if(LDType==1){
    out<<dW_ldpar<<"  "<<gamma_ldpar<<"  "<<ainf_ldpar<<"  "<<Delta_ldpar<<std::endl;
  }
  else if(LDType==2){
    out<<dW_ldpar<<"  "<<gamma_ldpar<<"  "<<ainf_ldpar<<"  "<<Delta_ldpar<<"  "<<Ex_ldpar<<"  "<<E0_ldpar<<"  "<<T_ldpar<<std::endl;
  }
  else if(LDType==3){
    out<<ainf_ldpar<<"  "<<Delta_ldpar<<std::endl;
  }
  out.precision(oldprc);
  out<<std::endl;
  
}


G4double G4NuDEXLevelDensity::GetNucleusTemperature(G4double ExcEnergy){

  if(!HasData){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(ExcEnergy<Ex_ldpar && LDType==2){
    return T_ldpar;
  }

  G4double Uval=ExcEnergy-Delta_ldpar;
  if(Uval<=0){return 0;}
  G4double a_ldpar=ainf_ldpar*(1.+dW_ldpar/Uval*(1.-std::exp(-gamma_ldpar*Uval)));
  if(LDType==3){
    a_ldpar=ainf_ldpar;
  }
  return std::sqrt(Uval/a_ldpar);


}


//Gilbert-Cameron:
G4double G4NuDEXLevelDensity::GetLevelDensity(G4double ExcEnergy,G4double spin,G4bool ,G4bool TotalLevelDensity){

  if(!HasData){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  //If A_Int even/odd --> spinx2 (spin_val*2) should be even/odd
  if(((A_Int+(G4int)(spin*2+0.01))%2)!=0 && (TotalLevelDensity==false)){
    return 0;
  }

  G4double Uval=ExcEnergy-Delta_ldpar;
  if(Uval<0){Uval=1.e-6;}

  //----------------------------------------------------------------
  // Back shifted: von Egidy et al., NP A481 (1988) 189
  if(LDType==3){
    G4double sig2=0.0888*std::pow(A_mass,2./3.)*std::sqrt(ainf_ldpar*Uval);
    G4double sig=std::sqrt(sig2);
    G4double rho=0.05893*std::exp(2.*std::sqrt(ainf_ldpar*Uval))/sig/std::pow(ainf_ldpar,0.25)/std::pow(Uval,1.25);
    G4double xj2=(spin+0.5)*(spin+0.5);
    G4double fj=(2.*spin+1.)/2./sig2*std::exp(-xj2/2./sig2);
    return 0.5*fj*rho;
  }
  //----------------------------------------------------------------


  //--------------------------------------------------------------------------------
  //statistical factor from eq. 39 of RIPL-3 manual, and sigma2 from eqs. 57-60
  G4double Uval_Sn=Sn-Delta_ldpar;
  G4double a_ldpar=ainf_ldpar*(1.+dW_ldpar/Uval*(1.-std::exp(-gamma_ldpar*Uval)));
  G4double a_ldpar_Sn=ainf_ldpar*(1.+dW_ldpar/Uval_Sn*(1.-std::exp(-gamma_ldpar*Uval_Sn)));
  G4double sigma2_f=0.01389*std::pow(A_mass,5./3.)/ainf_ldpar*std::sqrt(a_ldpar*Uval);
  G4double sigma2_f_Sn=0.01389*std::pow(A_mass,5./3.)/ainf_ldpar*std::sqrt(a_ldpar_Sn*Uval);
  G4double sigma2_d=(0.83*std::pow(A_mass,0.26))*(0.83*std::pow(A_mass,0.26));

  G4double sigma2;
  if(ExcEnergy<=Ed){
    sigma2=sigma2_d;//if ExcEnergy<Ed
  }
  else if(ExcEnergy<=Sn){
    sigma2=sigma2_d+(ExcEnergy-Ed)/(Sn-Ed)*(sigma2_f_Sn-sigma2_d);
  }
  else{
    sigma2=sigma2_f;
  }
  G4double statfactor=1./2.*(2.*spin+1.)/(2.*sigma2)*std::exp(-(spin+1/2.)*(spin+1/2.)/2./sigma2);
  if(TotalLevelDensity==true){
    statfactor=1;
  }
  //--------------------------------------------------------------------------------

  //CT + BSFG: Gilbert & Cameron, Can.J.Phys. 43 (1965) 1446
  if(LDType==2 && ExcEnergy<Ex_ldpar){
    G4double totalrho=std::exp((ExcEnergy-E0_ldpar)/T_ldpar)/T_ldpar;
    return totalrho*statfactor;
  }

  //Else: BSFGM (LDType==1 or ExcEnergy>Ex_ldpar)
  G4double rhotot_f=1./std::sqrt(2.*sigma2)/12.*std::exp(2.*std::sqrt(a_ldpar*Uval))/std::pow(a_ldpar,1./4.)/std::pow(Uval,5./4.);
  G4double rhotot_0=std::exp(1.)*a_ldpar/12./std::sqrt(sigma2)*std::exp(a_ldpar*Uval);
  G4double totalrho=1./(1./rhotot_f+1./rhotot_0);

  return totalrho*statfactor;

}


G4double G4NuDEXLevelDensity::EstimateInverse(G4double LevDen_iMeV,G4double spin,G4bool parity){

  //We assume that rho is a monotonically increasing function

  G4double tolerance=0.001; //the result will have this relative tolerance. 0.01 means 1%

  G4double xmin=0;
  G4double xmax=1;
  while(GetLevelDensity(xmax,spin,parity)<LevDen_iMeV){
    xmax*=2;
  }

  while(xmin/xmax<1-tolerance){
    G4double x0=(xmin+xmax)/2.;
    if(GetLevelDensity(x0,spin,parity)<LevDen_iMeV){
      xmin=x0;
    }
    else{
      xmax=x0;
    }
  }

  return (xmin+xmax)/2.;

}


G4double G4NuDEXLevelDensity::Integrate(G4double Emin,G4double Emax,G4double spin,G4bool parity){

  G4int nb=1000;
  G4double Integral=0,x1,x2,y1,y2;
  for(G4int i=0;i<nb;i++){
    x1=Emin+(Emax-Emin)*i/(G4double)(nb-1.);
    x2=Emin+(Emax-Emin)*(i+1.)/(G4double)(nb-1.);
    y1=GetLevelDensity(x1,spin,parity);
    y2=GetLevelDensity(x2,spin,parity);
    Integral+=(y1+y2)/2.*(x2-x1);
  }

  return Integral;
}

void G4NuDEXLevelDensity::PrintParameters(std::ostream &out){

  out<<" Level density type: "<<LDType<<std::endl;
  if(LDType==1){ // Back-Shifted-Fermi-Gas model
    out<<" ainf = "<<ainf_ldpar<<"  gamma = "<<gamma_ldpar<<"  dW = "<<dW_ldpar<<"  Delta = "<<Delta_ldpar<<std::endl;
  }
  else{
    out<<" ainf = "<<ainf_ldpar<<"  gamma = "<<gamma_ldpar<<"  dW = "<<dW_ldpar<<"  Delta = "<<Delta_ldpar<<"  T = "<<T_ldpar<<"  E0 = "<<E0_ldpar<<"  Ex = "<<Ex_ldpar<<std::endl;
  }

}




