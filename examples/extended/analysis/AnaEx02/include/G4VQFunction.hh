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
// $Id: G4VQFunction.hh,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
 
#ifndef G4VQFunction_h
#define G4VQFunction_h 1

#include "globals.hh"
#include "G4ios.hh"
//#include "Randomize.hh"
#include <vector>

class G4VQFunction
{
 public:

  G4VQFunction(std::vector<G4double>* AI, std::vector<G4double>* XI,
               std::vector<G4double>* MI=0, std::vector<G4double>* MA=0);
  virtual ~G4VQFunction();
      
 public:
  
  G4double               Y(std::vector<G4double>* df=0, std::vector<G4double>* da=0);//F&DF
  G4double               GetMin(G4int i)                    {return (*AMI)[i];}
  G4double               GetMax(G4int i)                    {return (*AMA)[i];}
  G4double               GetA(G4int i  )                    {return (*A)[i];}
  G4double               GetX(G4int i  )                    {return (*X)[i];}
  void                   SetMin(G4int i, G4double mi)       {(*AMI)[i]=mi;}
  void                   SetMax(G4int i, G4double ma)       {(*AMA)[i]=ma;}
  void                   SetA(G4int i  , G4double a )       {(*A)[i]  =a;}
  void                   SetX(G4int i  , G4double x )       {(*X)[i]  =x;}
  void SetMin(std::vector<G4double>* mi)
                {AMI->clear(); for(unsigned i=0;i<mi->size();i++)AMI->push_back((*mi)[i]);}
  void SetMax(std::vector<G4double>* ma)
                {AMA->clear(); for(unsigned i=0;i<ma->size();i++)AMA->push_back((*ma)[i]);}
  void SetA(std::vector<G4double>*    a)
                      {A->clear(); for(unsigned i=0;i<a->size();i++)A->push_back((*a)[i]);}
  void SetX(std::vector<G4double>*    x)
                      {X->clear(); for(unsigned i=0;i<x->size();i++)X->push_back((*x)[i]);}
  std::vector<G4double>* GetMinP()                          {return AMI;}
  std::vector<G4double>* GetMaxP()                          {return AMA;}
  std::vector<G4double>* GetAP()                            {return A;}
  std::vector<G4double>* GetXP()                            {return X;}
     
 private:

  virtual G4double Funct() =0;           // Function definition (pure virtual)
  virtual void Deriv(std::vector<G4double>* DF, std::vector<G4double>* DA); // Derivitives

 protected:

  std::vector<G4double>* A;               // Values of parameters
  std::vector<G4double>* X;               // Values of arguments
  std::vector<G4double>* AMI;             // Min values of parameters
  std::vector<G4double>* AMA;             // Min values of parameters
};

#endif
