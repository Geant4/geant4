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
// $Id: AnaEx02Function.hh,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
 
#ifndef AnaEx02Function_h
#define AnaEx02Function_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VQFunction.hh"
//#include "Randomize.hh"
#include <vector>

class AnaEx02Function : public G4VQFunction
{
 public:

  AnaEx02Function(std::vector<G4double>* AI, std::vector<G4double>* XI,
               std::vector<G4double>* MI=0, std::vector<G4double>* MA=0);
  //~AnaEx02Function();
      
 public:
  
  //G4double GetA() {return a;}
     
 private:

  G4double Funct();           // Function definition for the particular application
  //void Deriv(std::vector<G4double>* DF, std::vector<G4double>* PL);// DerivitivesOverParams

};

#endif
