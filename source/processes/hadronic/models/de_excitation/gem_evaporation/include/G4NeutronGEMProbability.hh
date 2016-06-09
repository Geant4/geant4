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
// $Id: G4NeutronGEMProbability.hh,v 1.2 2004/12/07 13:46:49 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4NeutronGEMProbability_h
#define G4NeutronGEMProbability_h 1


#include "G4GEMProbability.hh"


class G4NeutronGEMProbability : public G4GEMProbability
{
public:
    // Only available constructor
    G4NeutronGEMProbability();
		
    ~G4NeutronGEMProbability() {}
private:  
    // Copy constructor
    G4NeutronGEMProbability(const G4NeutronGEMProbability &right);

    const G4NeutronGEMProbability & operator=(const G4NeutronGEMProbability &right);
    G4bool operator==(const G4NeutronGEMProbability &right) const;
    G4bool operator!=(const G4NeutronGEMProbability &right) const;
  

private:
    
    virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
        {
            return 0.76+2.2/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),1.0/3.0);
        }
	
    virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
        {
            return (2.12/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),2.0/3.0)-0.05)*MeV/
	      CalcAlphaParam(fragment);
        }
    
    // Excitation energy levels 
    std::vector<G4double> ExcitEnergies;
    // Spin of excitation energy levels 
    std::vector<G4double> ExcitSpins;

    std::vector<G4double> ExcitLifetimes;
    
};


#endif
