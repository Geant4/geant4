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
// $Id: NeutronEvapProb.hh,v 1.1 2003-10-08 12:32:18 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4NeutronEvaporationProbability_h_
#define G4NeutronEvaporationProbability_h_ 1


#include "MyEvapProb.hh"


class NeutronEvapProb : public MyEvapProb
{
public:
  // Only available constructor
  NeutronEvapProb();
		
  ~NeutronEvapProb() {}
  virtual G4double CoulombBarrier(G4Fragment& fragm){
    return 0;
  }
private:  
  // Copy constructor
  NeutronEvapProb(const NeutronEvapProb &right);

  const NeutronEvapProb & operator=(const NeutronEvapProb &right);
  G4bool operator==(const NeutronEvapProb &right) const;
  G4bool operator!=(const NeutronEvapProb &right) const;
  

private:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
  { return 0.76+2.2/pow(G4double(fragment.GetA()-GetA()),1.0/3.0);}
	
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
  { return (2.12/pow(G4double(fragment.GetA()-GetA()),2.0/3.0) - 0.05)*MeV/CalcAlphaParam(fragment); }

  // Excitation energy levels 
  std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  std::vector<G4int> ExcitSpins;

};


#endif
