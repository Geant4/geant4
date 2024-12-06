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



#include "G4NuDEXInternalConversion.hh"



//If alpha>0, use that value
G4bool G4NuDEXInternalConversion::SampleInternalConversion(G4double Ene,G4int multipolarity,G4double alpha,G4bool CalculateProducts){

  if(theZ<MINZINTABLES){ //then we have no info
    if(alpha<0){
      Ne=0;
      Ng=0;
      return false;
    }
    else{
      G4double rand=theRandom4->Uniform(0,alpha+1);
      if(rand<alpha){ //then electron conversion
	Ne=1;
	Ng=0;
	Eele[0]=Ene; //which is not correct, but we don't know the binding energy
	return true;
      }
      return false;
    }
  }


  Ne=0;
  Ng=0;

  if(multipolarity==0){ //maybe it is better to return true ... ?? --> no
    //return true;
    if(alpha<=0){
      return false;
    }
    //NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  G4bool usegivenalpha=true;
  if(NShells==0 || std::abs(multipolarity)>ICC_NMULTIP){return false;}
  if(alpha<0){
    usegivenalpha=false;
    alpha=GetICC(Ene,multipolarity);
  }

  G4double rand=theRandom4->Uniform(0,alpha+1);
  if(rand<alpha){ //then electron conversion
    if(!CalculateProducts){return true;}
    //Select the orbital:
    if(usegivenalpha){rand=rand*GetICC(Ene,multipolarity)/alpha;} //renormalize rand to our alpha
    G4double cumul=0;
    for(G4int i=1;i<NShells;i++){
      cumul+=GetICC(Ene,multipolarity,i);
      //std::cout<<Ene<<"  "<<multipolarity<<"  "<<i<<"  "<<GetICC(Ene,multipolarity,i)<<"  "<<rand-1<<std::endl;
      if(cumul>=rand || multipolarity==0){ //then is this orbital
	Ne=1;
	Eele[0]=Ene-BindingEnergy[i];
	FillElectronHole(i); //now there is a hole there, in the filling procedure we emitt gammas and/or electrons
	if(Eele[0]<0){
	  std::cout<<" For Z = "<<theZ<<" and orbital "<<OrbitalName[i]<<" --> Ene = "<<Ene<<" and BindingEnergy = "<<BindingEnergy[i]<<std::endl;
	  std::cout<<" Given alpha is "<<alpha<<" ("<<usegivenalpha<<"), rand = "<<rand<<" and tabulated alpha for Ene = "<<Ene<<" and mult = "<<multipolarity<<" is "<<GetICC(Ene,multipolarity)<<" -- cumul = "<<cumul<<std::endl;
	  for(G4int j=1;j<=NShells;j++){
	    std::cout<<j<<"  "<<GetICC(Ene,multipolarity,j)<<std::endl;
	  }
	  Eele[0]=0;
	}
	return true;
      }
    }
    std::cout<<" ############ Warning in "<<__FILE__<<", line "<<__LINE__<<" ############"<<std::endl;
    std::cout<<" Given alpha is "<<alpha<<" and tabulated alpha for Ene = "<<Ene<<" and mult = "<<multipolarity<<" is "<<GetICC(Ene,multipolarity)<<" -- cumul = "<<cumul<<std::endl;
    for(G4int i=1;i<=NShells;i++){
      std::cout<<i<<"  "<<GetICC(Ene,multipolarity,i)<<std::endl;
    }
    Ne=1;
    Eele[0]=Ene-BindingEnergy[NShells-1];
    return true;
  }

  return false;
}


void G4NuDEXInternalConversion::FillElectronHole(G4int i_shell){

  //A very simplified version of the process (... and false). It can be done with accuracy with G4AtomicTransitionManager

  G4double fluoyield=0;
  if(i_shell==1){ //K-shell
    //Hubbell et al. (1994) formula for the fluorescence yield:
    G4double C0=0.0370,C1=0.03112,C2=5.44e-5,C3=-1.25e-6;
    G4double w_fac=std::pow(C0+C1*theZ+C2*theZ*theZ+C3*theZ*theZ*theZ,4);
    fluoyield=w_fac/(1.+w_fac);
  }
  else if(i_shell>=2 && i_shell<=4){ //L-shell
    //Hubbell et al. (1994) formula for the fluorescence yield:
    if(theZ>=3 && theZ<=36){
      fluoyield=1.939e-8*std::pow(theZ,3.8874);
    }
    else if(theZ>36){
      G4double C0=0.17765,C1=0.00298937,C2=8.91297e-5,C3=-2.67184e-7;
      G4double w_fac=std::pow(C0+C1*theZ+C2*theZ*theZ+C3*theZ*theZ*theZ,4);
      fluoyield=w_fac/(1.+w_fac);
    }
  }


  G4double rand=theRandom4->Uniform(0,1);
  if(rand<fluoyield){ //gamma emission
    Egam[Ng]=BindingEnergy[i_shell];
    Ng++;
  }
  else{ //electron emission
    Eele[Ne]=BindingEnergy[i_shell];
    Ne++;
  }


}







//If i_shell<0 --> the total alpha
G4double G4NuDEXInternalConversion::GetICC(G4double Ene,G4int multipolarity,G4int i_shell){

  if(theZ<MINZINTABLES){ //then we have no info
    return 0;
  }

  if(NShells==0 || std::abs(multipolarity)>ICC_NMULTIP){return 0;}

  //-----------------------------------------
  //Total:
  //The following line does not work, due to interpolation below binding energies:
  //if(i_shell<0){i_shell=NShells;}

  if(i_shell<0){
    G4double result=0;
    for(G4int i=1;i<NShells;i++){
      result+=GetICC(Ene,multipolarity,i);
    }
    return result;
  }
  //-----------------------------------------


  if(Ene<BindingEnergy[i_shell]){return 0;}

  if(np[i_shell]==0){
    std::cout<<" shell "<<i_shell<<" has not been initialized"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  if(i_shell==NShells && Ene<Eg[i_shell][0]){ //then we cannot extrapolate, because of the binding energies of the different shells
    G4double total=0;
    for(G4int i=1;i<NShells;i++){total+=GetICC(Ene,multipolarity,i);}
    return total;
  }

  if(multipolarity>0){
    return Interpolate(Ene,np[i_shell],Eg[i_shell],Icc_E[multipolarity-1][i_shell]);
  }
  else if(multipolarity<0){
    return Interpolate(Ene,np[i_shell],Eg[i_shell],Icc_M[(-multipolarity)-1][i_shell]);
  }
  return 0;
}


G4NuDEXInternalConversion::G4NuDEXInternalConversion(G4int Z){
  theZ=Z;
  NShells=0;
  Ne=Ng=0;
  for(G4int i=0;i<ICC_MAXNSHELLS;i++){
    Eg[i]=0; np[i]=0;  BindingEnergy[i]=0;
    for(G4int j=0;j<ICC_NMULTIP;j++){
      Icc_E[j][i]=0; Icc_M[j][i]=0; 
    }
  }
  theRandom4= new G4NuDEXRandom(1234567);
}


G4NuDEXInternalConversion::~G4NuDEXInternalConversion(){
  for(G4int i=0;i<ICC_MAXNSHELLS;i++){
    if(Eg[i]!=0){delete [] Eg[i];}
    for(G4int j=0;j<ICC_NMULTIP;j++){
      if(Icc_E[j][i]!=0){delete [] Icc_E[j][i];}
      if(Icc_M[j][i]!=0){delete [] Icc_M[j][i];}
    }
  }
  delete theRandom4;
}

void G4NuDEXInternalConversion::PrintICC(std::ostream &out){

  char word[1000];
  out<<" ######################################################################################################################################### "<<std::endl;
  out<<" ICC"<<std::endl;
  out<<" Z = "<<theZ<<std::endl;
  out<<" NShells = "<<NShells<<std::endl;
  out<<" ----------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  out<<" Total calculated from the sum of the partials:"<<std::endl;
  out<<"     E_g         E1          E2          E3          E4          E5          M1          M2          M3          M4          M5 "<<std::endl;
  for(G4int j=0;j<np[NShells];j++){
    snprintf(word,1000,"%10.4g",Eg[NShells][j]); out<<word;
    for(G4int k=0;k<ICC_NMULTIP;k++){
      snprintf(word,1000,"  %10.4g",Icc_E[k][NShells][j]); out<<word;
    }
    for(G4int k=0;k<ICC_NMULTIP;k++){
      snprintf(word,1000,"  %10.4g",Icc_M[k][NShells][j]); out<<word;
    }
    out<<std::endl;
  }
  out<<" ----------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  for(G4int i=0;i<NShells;i++){
    out<<" ----------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    out<<" Binding energy = "<<BindingEnergy[i]<<" MeV -  OrbitalName = "<<OrbitalName[i]<<" -  np = "<<np[i]<<std::endl;
    out<<"     E_g         E1          E2          E3          E4          E5          M1          M2          M3          M4          M5 "<<std::endl;
    for(G4int j=0;j<np[i];j++){
      snprintf(word,1000,"%10.4g",Eg[i][j]); out<<word;
      for(G4int k=0;k<ICC_NMULTIP;k++){
	snprintf(word,1000,"  %10.4g",Icc_E[k][i][j]); out<<word;
      }
      for(G4int k=0;k<ICC_NMULTIP;k++){
	snprintf(word,1000,"  %10.4g",Icc_M[k][i][j]); out<<word;
      }
      out<<std::endl;
    }
    out<<" ----------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  }
  out<<" ########################################################################################################################################## "<<std::endl;
}

void G4NuDEXInternalConversion::Init(const char* fname){

  if(theZ<MINZINTABLES){ //then we have no info
    return;
  }

  if(NShells!=0){ //Init only once
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  std::ifstream in(fname);
  if(!in.good()){
    std::cout<<" ################ Error opening "<<fname<<" ################"<<std::endl;
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  std::string word;
  NShells=1;
  while(in>>word){
    if(word.c_str()[0]=='Z' && word.c_str()[1]=='='){
      if(std::atoi(&(word.c_str()[2]))==theZ){
	in>>word>>word;
	G4int orbindex;
	if(word==std::string("Total")){
	  in.ignore(1000,'\n');
	  in.ignore(1000,'\n');
	  orbindex=0;
	}
	else{
	  orbindex=NShells;
	  in>>word>>word>>BindingEnergy[NShells]>>word;
	  BindingEnergy[NShells]*=1.e-6; // all in MeV
	  in.ignore(1000,'\n');
	  in.ignore(1000,'\n');
	}
	//--------------------------------------------------------------------------------
	size_t sz,sz2;
	G4int np_tmp=0;
	G4double Eg_tmp[1000],Icc_E_tmp[ICC_NMULTIP][100],Icc_M_tmp[ICC_NMULTIP][100];
	while(getline(in,word)){
	  if(word.size()<100){
	    np[orbindex]=np_tmp;
	    Eg[orbindex]=new G4double[np_tmp];
	    for(G4int j=0;j<np_tmp;j++){
	      Eg[orbindex][j]=Eg_tmp[j];
	    }
	    for(G4int i=0;i<ICC_NMULTIP;i++){
	      Icc_E[i][orbindex]=new G4double[np_tmp];
	      Icc_M[i][orbindex]=new G4double[np_tmp];
	    }
	    for(G4int i=0;i<ICC_NMULTIP;i++){
	      for(G4int j=0;j<np_tmp;j++){
		Icc_E[i][orbindex][j]=Icc_E_tmp[i][j];
		Icc_M[i][orbindex][j]=Icc_M_tmp[i][j];
	      }
	    }
	    if(orbindex!=0){NShells++;}
	    break;
	  }
	  else{
	    sz=0;
	    Eg_tmp[np_tmp]=std::stof(word,&sz2); 
	    Eg_tmp[np_tmp]*=1.e-3; //all to MeV
	    sz+=sz2;
	    for(G4int i=0;i<ICC_NMULTIP;i++){
	      Icc_E_tmp[i][np_tmp]=std::stof(word.substr(sz),&sz2); sz+=sz2; 
	    }
	    for(G4int i=0;i<ICC_NMULTIP;i++){
	      Icc_M_tmp[i][np_tmp]=std::stof(word.substr(sz),&sz2); sz+=sz2; 
	    }
	    if((G4int)(std::stof(word.substr(sz),&sz2)+0.01)!=theZ){ 
	      NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
	    }
	    sz+=sz2;
	    sz2=word.find_first_not_of(' ',sz);
	    if(np_tmp==0){OrbitalName[orbindex]=word.substr(sz2,word.substr(sz2).size()-1);}
	    np_tmp++;
	  }
	}
	if(orbindex==0){break;}
	//--------------------------------------------------------------------------------
      }
    }
  }
  if(!in.good()){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }
  in.close();

  MakeTotal();
}



// Total Icc goes to nShell=NShells
void G4NuDEXInternalConversion::MakeTotal(){

  if(np[0]==0 || Eg[0]==0){
    NuDEXException(__FILE__,std::to_string(__LINE__).c_str(),"##### Error in NuDEX #####");
  }

  //We evaluate it in the same frame as the total given by the data:
  BindingEnergy[NShells]=0;
  np[NShells]=np[0];
  Eg[NShells]=new G4double[np[NShells]];
  for(G4int i=0;i<ICC_NMULTIP;i++){
    Icc_E[i][NShells]=new G4double[np[NShells]];
    Icc_M[i][NShells]=new G4double[np[NShells]];
  }
  for(G4int k=0;k<np[NShells];k++){
    for(G4int j=0;j<ICC_NMULTIP;j++){
      Icc_E[j][NShells][k]=0;
      Icc_M[j][NShells][k]=0;
    }
  }
  

  for(G4int k=0;k<np[NShells];k++){
    Eg[NShells][k]=Eg[0][k];
    for(G4int i=1;i<NShells;i++){
      for(G4int j=0;j<ICC_NMULTIP;j++){
	Icc_E[j][NShells][k]+=GetICC(Eg[NShells][k],j+1,i);
	Icc_M[j][NShells][k]+=GetICC(Eg[NShells][k],-j-1,i);
      }
    }
  }
  
}

//if val>x[npmax] then ---> return 0
G4double G4NuDEXInternalConversion::Interpolate(G4double val,G4int npoints,G4double* x,G4double* y){

  G4int i_interp=-1;
  for(G4int i=1;i<npoints;i++){
    if(x[i]>=val){i_interp=i-1; break;}
  }
  if(i_interp<0){return 0;}


  /*
  //y=a0+a1*x
  G4double a1=(y[i_interp+1]-y[i_interp])/(x[i_interp+1]-x[i_interp]);
  G4double a0=y[i_interp]-a1*x[i_interp];

  return (a0+a1*val);
  */

  //It is better a log-log interpolation:
  if(y[i_interp+1]<=0 || y[i_interp]<=0 || x[i_interp+1]<=0 || x[i_interp]<=0){
    //y=a0+a1*x
    G4double a1=(y[i_interp+1]-y[i_interp])/(x[i_interp+1]-x[i_interp]);
    G4double a0=y[i_interp]-a1*x[i_interp];

    return (a0+a1*val);
  }

  //log(y)=a0+a1*log(x)
  G4double a1=std::log(y[i_interp+1]/y[i_interp])/std::log(x[i_interp+1]/x[i_interp]);
  G4double a0=std::log(y[i_interp])-a1*std::log(x[i_interp]);

  G4double result=std::exp(a0+a1*std::log(val));
  return result;
}
