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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: MyEvapProb.hh,v 1.1 2003-10-08 12:32:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4EvaporationProbability_h_
#define G4EvaporationProbability_h_ 1


#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "MyLDP.hh"


class MyEvapProb : public G4VEmissionProbability
{
public:
  // Only available constructor
  MyEvapProb(const G4int anA, const G4int aZ, const G4double aGamma);

  ~MyEvapProb() 
  {
    if (theEvapLDPptr != 0) delete theEvapLDPptr;
  }


	
  G4double GetZ(void) const { return theZ; }
	
  G4double GetA(void) const { return theA;} 

protected:
  G4double BindingEnergy();

  void SetExcitationEnergiesPtr(std::vector<G4double> * anExcitationEnergiesPtr) 
  {ExcitationEnergies = anExcitationEnergiesPtr;}

  void SetExcitationSpinsPtr(std::vector<G4int> * anExcitationSpinsPtr)
  {ExcitationSpins = anExcitationSpinsPtr;}

  
  // Default constructor
  MyEvapProb():theEvapLDPptr(NULL) {}
private:
  // Copy constructor
  MyEvapProb(const MyEvapProb &right);

  const MyEvapProb & operator=(const MyEvapProb &right);
  G4bool operator==(const MyEvapProb &right) const;
  G4bool operator!=(const MyEvapProb &right) const;
  
public:
  G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy);
  G4double SampleEnergy(const G4Fragment& fragment,const G4double anEnergy);

private:

  G4double CalcProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy);
  virtual G4double CCoeficient(const G4double aZ) const {return 0.0;};

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const {return 1.0;}
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const {return 1.0;}
  virtual G4double CoulombBarrier(const G4Fragment& fragm){return 0;};
  G4double FindMaximum(G4double left,G4double right,G4double hint=0);

  double ro(double U,int Prime=0,int mat=0);
  double func(double U,int Prime=0);
  double Integrator(double Min,double Max,double Left,double Right);
  

  // Data Members

  G4VLevelDensityParameter * theEvapLDPptr;
	
  G4int theA;
  G4int theZ;
  G4Fragment fragm;
  G4double m_NuclA;
  G4double m_NuclB;
  G4double Beta;
  G4double Emax;
  G4double Estar;
  G4double m_delta0;

  // Gamma is A_f(2S_f+1) factor, where A_f is fragment atomic 
  // number and S_f is fragment spin
  G4double Gamma;

  // Discrete Excitation Energies 
  std::vector<G4double> * ExcitationEnergies;

  //
  std::vector<G4int> * ExcitationSpins;

};


#endif
