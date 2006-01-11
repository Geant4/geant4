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
// $Id: G4QFunFit.hh,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
 
#ifndef G4QFunFit_h
#define G4QFunFit_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ExceptionHandler.hh"
#include "G4VQFunction.hh"
//#include "Randomize.hh"
#include <vector>
#include "G4Timer.hh"
#include "time.h"

class G4QFunFit
{
 public:

  G4QFunFit(G4bool fl, std::vector<G4double>* AI, std::vector<G4double>* DI,
          		std::vector<G4double>* EX, G4VQFunction* usefun, G4int Np,
            std::vector<G4double>* MI=0, std::vector<G4double>* MA=0);

  ~G4QFunFit();
      
 public:
  // S=Hi2/2(out), N1=Red(in), N2=Inc(in), N3=maxIter#(in), EPS=accuracy(in), IT=output(in)
  G4double FIT(G4int N1, G4int N2, G4int N3, G4double EPS, G4int IT); // ifError -> S=-S
  //
  G4double SGZ();             // Return S, but fill Z0 as well
  void MCONV(G4int N);
  void ERRORF();              // Additional (more accurate) calculation of errors
  std::vector<G4double>* GetAP()     {return A;}
  std::vector<G4double>* GetDP()     {return DA;}
  std::vector<G4double>* GetSIG()    {return SIGMA;}
     
 private:

  void MEXIT(G4int IR, G4int N, G4int II); // additional slave function
  G4double SCAL();                   // Calculates ER = variance of theoretical function Y
  void MONITO(G4double S, G4int NN3, G4int IT, G4int ENDFLG,
              G4double GT,G4double AKAPPA,G4double ALAMBD);

 private:
  G4int NA;                    // a#of parameters
  G4int NS;                    // a#of experimental point (@@ take from data)
  G4int NP;                    // Length of the Y,DY,X() line (@@ take from data)
  G4int fixedpr;               // fixed parameter number INDFLG(1)-1
  G4bool fulik;                // fulik=false -> hi2, fulik=true -> Liklyho
  //G4bool hasafix;              // a flag of existing fixed parameterod
  G4double AM;                 // Primary big number
  G4double RP;                 // Primary small number
  G4double AMM;                // Primary intermediate
  G4double APR;                // squared big number
  G4double APS;                // big number
  G4double AP;                 // small number
  std::vector<G4double> PL;    // Length Lpar
  std::vector<G4double> Z;     // Length L3par
  std::vector<G4double> G;     // Length Lpar
  std::vector<G4double>* DF;   // Length Lpar
  std::vector<G4double>* A;    // Length Lpar
  std::vector<G4double>* DA;   // Proposed change of parameters
  std::vector<G4double>* PL0;  // Length Lpar INPUT start values for search ranges
  std::vector<G4double>* AMN;  // Length Lpar
  std::vector<G4double>* AMX;  // Length Lpar
  std::vector<G4double>* SIGMA;// Length Lpar
  std::vector<G4double>* EXDA; // Length LLpar
  std::vector<G4double> DO;    // Copy of the old Proposed change of parameters
  std::vector<G4double> R;     // Length Lpar
  std::vector<G4double> ERROR; // Length L1par
  std::vector<G4double> Z0;    // Length L3par
  std::vector<G4int>    JR;    // PL-reduced array of pointers to the not PL-reduced (full)
  G4VQFunction* fumfun;        // Pointer to Y-DY function  
  G4bool miF;                  // default AMI creation
  G4bool maF;                  // default AMA creation
};

#endif
