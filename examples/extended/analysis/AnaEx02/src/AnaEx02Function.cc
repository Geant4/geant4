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
// $Id: AnaEx02Function.cc,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------- 

//#define debug
//#define edebug
//#define vdebug

#include "AnaEx02Function.hh"

AnaEx02Function::AnaEx02Function(std::vector<G4double>* AI, std::vector<G4double>* XI,
                           std::vector<G4double>* MI, std::vector<G4double>* MA)
                               : G4VQFunction(AI,XI,MI,MA) {}

G4double AnaEx02Function::Funct()
{
#ifdef debug
  G4cout<<"AnaEx02Function::Funct: X0="<<(*X)[0]<<",X1="<<(*X)[1]<<",X2="<<(*X)[2]<<",X3="
        <<(*X)[3]<<",A0="<<(*A)[0]<<",A1="<<(*A)[1]<<",A2="<<(*A)[2]<<",A3="<<(*A)[3]
        <<",A4="<<(*A)[4]<<",A5="<<(*A)[5]<<G4endl;
#endif
  return (*A)[5]*(*X)[3]+(*A)[4]*std::exp((*A)[3]*(*X)[2])*(std::exp((*A)[2]*(*X)[1])+
 (*A)[1]*std::log((*A)[0]*(*X)[0]));
}

// A default 1%-of-var(PL) calculation of derivitives
//void AnaEx02Function::Deriv(std::vector<G4double>* DF, std::vector<G4double>* PL)
//{
//#ifdef debug
//  G4cout<<"AnaEx02Function::Deriv: X0="<<(*X)[0]<<",X1="<<(*X)[1]<<",X2="<<(*X)[2]<<",X3="
//        <<(*X)[3]<<",A0="<<(*A)[0]<<",A1="<<(*A)[1]<<",A2="<<(*A)[2]<<",A3="<<(*A)[3]
//        <<",A4="<<(*A)[4]<<",A5="<<(*A)[5]<<G4endl;
//#endif
//  G4double a2=std::exp((*A)[2]*(*X)[1]);
//  G4double a3=std::exp((*A)[3]*(*X)[2]);
//  G4double al=std::log((*A)[5]*(*X)[3]);
//		G4double sl=a3+(*A)[4]*al;
//  G4double s2=a2*sl;
//  G4double aa=(*A)[1]*a2;
//  DF->clear();
//  if((*PL)[0]>0) DF->push_back((*X)[0]);
//		else DF->push_back(0.);
//  if((*PL)[1]>0) DF->push_back(s2);
//		else DF->push_back(0.);
//  if((*PL)[2]>0) DF->push_back((*X)[0]*(*A)[1]*s2);
//		else DF->push_back(0.);
//  if((*PL)[3]>0) DF->push_back((*X)[1]*aa*a3);
//		else DF->push_back(0.);
//  if((*PL)[4]>0) DF->push_back(aa*al);
//		else DF->push_back(0.);
//  if((*PL)[5]>0) DF->push_back(aa*(*A)[4]/(*A)[5]);
//		else DF->push_back(0.);
//#ifdef debug
//  G4cout<<"AnaEx02Function::Deriv: done,D0="<<(*DF)[0]<<",D1="<<(*DF)[1]<<",D2="<<(*DF)[2]
//        <<",D3="<<(*DF)[3]<<",D4="<<(*DF)[4]<<",D5="<<(*DF)[5]<<G4endl;
//#endif
//  return;
//}
