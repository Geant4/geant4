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
//
// $Id: G4VQFunction.cc,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------- 

//#define debug
//#define ddebug
#define edebug
//#define vdebug

#include "G4VQFunction.hh"

G4VQFunction::G4VQFunction(std::vector<G4double>* AI, std::vector<G4double>* XI,
                           std::vector<G4double>* MI, std::vector<G4double>* MA)
{
  static const G4double BIGV=1.e227;          // Big value
  G4int i;                                    // Prototype of the index
  G4int NX=XI->size();                        // #of arguments
  X = new std::vector<G4double>;              // create a new vector of arguments
  for(i=0; i<NX; i++) X->push_back((*XI)[i]); // fill arguments
  G4int NA=AI->size();                        // #of parameters
  A = new std::vector<G4double>;              // create a new vector of parameters
  for(i=0; i<NA; i++) A->push_back((*AI)[i]); // fill parameters
#ifdef debug
  G4cout<<"G4QFunction::: Created NA="<<NA<<", NX="<<NX<<G4endl;
#endif
  AMI = new std::vector<G4double>;            // create a new vector of minimum values
  if(!MI) for(i=0; i<NA; i++) AMI->push_back(-BIGV); // fill minimum by the default values
  else
		{
    G4int NMI=MI->size();                     // #of parameters
#ifdef edebug
    if(NMI!=NA) G4cout<<"%Warning% G4QFunction::: NA="<<NA<<" # NMIN="<<NMI<<G4endl;
#endif
    for(i=0; i<NMI; i++) AMI->push_back((*MI)[i]); // fill minimum
  }
  AMA = new std::vector<G4double>;            // create a new vector of minimum values
  if(!MA) for(i=0; i<NA; i++) AMA->push_back(BIGV); // fill maximum by the default values
  else
		{
    G4int NMA=MA->size();                     // #of parameters
#ifdef edebug
    if(NMA!=NA) G4cout<<"%Warning% G4QFunction::: NA="<<NA<<" # NMAX="<<NMA<<G4endl;
#endif
    for(i=0; i<NMA; i++) AMA->push_back((*MA)[i]); // fill maximum
  }
}

G4VQFunction::~G4VQFunction()                 // Always clean up what was created
{
  delete X;
  delete A;
  delete AMI;
  delete AMA;
}

G4double G4VQFunction::Y(std::vector<G4double>* df, std::vector<G4double>* da) // da=0 -> Y
{
#ifdef debug
  G4cout<<"G4QFunction::Y: is called, da="<<da<<", df="<<df<<G4endl; 
#endif
  if(da) Deriv(df, da);
  G4double y= Funct();
#ifdef debug
  G4cout<<"G4VQFunction::Y: done, Y="<<y;
  if(da) G4cout<<",D0="<<(*df)[0]<<",D1="<<(*df)[1]<<",D2="<<(*df)[2]<<",D3="<<(*df)[3]
               <<",D4="<<(*df)[4]<<",D5="<<(*df)[5];
  G4cout<<G4endl;
#endif
  return y;
}

// A default 1%-of-var(PL) calculation of derivitives
void G4VQFunction::Deriv(std::vector<G4double>* DF, std::vector<G4double>* PL)
{
  static const G4double RP=1.e-6;                        // small number (@@)
  G4int NA=A->size();                                    // #of parameters
  DF->clear();                                           // Empty the existing values
#ifdef ddebug
  G4int NX=X->size();                                    // #of arguments
  G4cout<<"G4QFunction::Deriv: Called, NX="<<NX<<", NA="<<NA<<G4endl; 
  for(G4int J=0; J<NA; J++) G4cout<<",A["<<J<<"]="<<(*A)[J]; 
  for(G4int J=0; J<NX; J++) G4cout<<",X["<<J<<"]="<<(*X)[J]; 
  G4cout<<G4endl;
#endif
  G4double Y=Funct();                                    // Function Value
#ifdef debug
  G4cout<<"G4QFunction::Deriv: Y="<<Y<<", NA="<<NA<<G4endl; 
#endif
  for(G4int I=0; I<NA; I++)
		{
    if((*PL)[I]<=0) DF->push_back(0.);
    else
    {
      G4double AI=(*A)[I];
      G4double HI=0.01*(*PL)[I];                         // 1%PL Step (dx_i)
      G4double DI=RP*std::fabs(AI);                      // RP*A Step (dx_i)
      if(HI<DI) HI=DI;                                   // Select the smallest
      G4double AJ=AI+HI;                                 // Second parameter value
						if(AJ>(*AMA)[I])                                   // ... is more than the Max Limit
						{
        AJ=AI-HI;
        HI=-HI;
        if(AJ<(*AMI)[I])                                 // ... is less than the Min Limit
								{
          G4double HL=AI-(*AMI)[I];                      // lower limit
          G4double HU=(*AMA)[I]-AI;                      // upper limit
          if(HL>HU) HI=-HL/2;                            // more in lower part
          else      HI=HU/2;                             // more in upper part
          AJ=AI+HI;
        }
      }
      (*A)[I]=AJ;
#ifdef debug
      G4cout<<"G4QFunction::Deriv: par # "<<I; 
      for(G4int J=0; J<NX; J++) G4cout<<",X["<<J<<"]="<<(*X)[J]; 
      G4cout<<G4endl;
#endif
      //#ifdef edebug
      if(!HI) G4cout<<"%Warning% G4QFunction::Deriv: HI=0"<<G4endl; 
      //#endif
      G4double YJ=Funct();
      G4double Da=(YJ-Y)/HI;
      DF->push_back(Da);
      (*A)[I]=AI;                                        // Recover the parameter value
#ifdef ddebug
      G4cout<<"G4QFunction::Deriv: D["<<I<<"]="<<Da<<",A="<<(*A)[I]<<",B="<<AJ<<",Y="<<Y
            <<",Z="<<YJ<<",dY="<<YJ-Y<<",dA="<<HI<<"="<<AJ-(*A)[I]<<",X="<<(*X)[I]<<G4endl;
#endif
    }
  }
#ifdef debug
  G4cout<<"G4QFunction::Deriv: *Done* NA="<<NA; 
  for(G4int J=0; J<NA; J++) G4cout<<",D["<<J<<"]="<<(*DF)[J]; 
  G4cout<<G4endl;
#endif
  return;
}
