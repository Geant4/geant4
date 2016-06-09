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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001) 
//



#ifndef G4GEMProbability_h
#define G4GEMProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4VCoulombBarrier.hh"
#include "G4PairingCorrection.hh"

class G4GEMProbability : public G4VEmissionProbability
{
public:
    // Only available constructor
    G4GEMProbability(const G4int anA, const G4int aZ, const G4double aSpin) : 
        theA(anA), theZ(aZ), Spin(aSpin), theCoulombBarrierPtr(0),
        ExcitationEnergies(0), ExcitationSpins(0), ExcitationLifetimes(0), Normalization(1.0)
        {
            theEvapLDPptr = new G4EvaporationLevelDensityParameter;
        }
    
    ~G4GEMProbability()
        {
            if (theEvapLDPptr != 0) delete theEvapLDPptr;
        }


	
    G4double GetZ(void) const { return theZ; }
	
    G4double GetA(void) const { return theA;}

    G4double GetSpin(void) const { return Spin; }

    G4double GetNormalization(void) const { return Normalization; }
    
    void SetCoulomBarrier(const G4VCoulombBarrier * aCoulombBarrierStrategy)
        {
            theCoulombBarrierPtr = aCoulombBarrierStrategy;
        }

    G4double GetCoulombBarrier(const G4Fragment& fragment) const 
        {
            if (theCoulombBarrierPtr) 
	      {
		G4int Acompound = static_cast<G4int>(fragment.GetA());
		G4int Zcompound = static_cast<G4int>(fragment.GetZ());
		return theCoulombBarrierPtr->GetCoulombBarrier(Acompound-theA, Zcompound-theZ,
							       fragment.GetExcitationEnergy()-
							       G4PairingCorrection::GetInstance()->
							       GetPairingCorrection(Acompound,Zcompound));
	      }
            else 
	      {
		return 0.0;
	      }
        }
    

    virtual G4double CalcAlphaParam(const G4Fragment & ) const {return 1.0;}
    virtual G4double CalcBetaParam(const G4Fragment & ) const {return 1.0;}
    
protected:
  
    void SetExcitationEnergiesPtr(std::vector<G4double> * anExcitationEnergiesPtr) 
        {
            ExcitationEnergies = anExcitationEnergiesPtr;
        }
  
    void SetExcitationSpinsPtr(std::vector<G4double> * anExcitationSpinsPtr)
        {
            ExcitationSpins = anExcitationSpinsPtr;
        }

    void SetExcitationLifetimesPtr(std::vector<G4double> * anExcitationLifetimesPtr)
        {
            ExcitationLifetimes = anExcitationLifetimesPtr;
        }

    void SetCoulombBarrierStrategy(G4VCoulombBarrier * aCoulombBarrier)
        {
            theCoulombBarrierPtr = aCoulombBarrier;
        }

  
    // Default constructor
    G4GEMProbability() {}
private:
    // Copy constructor
    G4GEMProbability(const G4GEMProbability &right);
    
    const G4GEMProbability & operator=(const G4GEMProbability &right);
    G4bool operator==(const G4GEMProbability &right) const;
    G4bool operator!=(const G4GEMProbability &right) const;
    
public:
    G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy);
  
private:

    G4double CalcProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy,
                             const G4double V);
    virtual G4double CCoeficient(const G4double ) const {return 0.0;};

    
    G4double I0(const G4double t);
    G4double I1(const G4double t, const G4double tx);
    G4double I2(const G4double s, const G4double sx);
    G4double I3(const G4double s, const G4double sx);
    
    // Data Members
    
    G4VLevelDensityParameter * theEvapLDPptr;
	
    G4int theA;
    G4int theZ;
    
    // Spin is fragment spin
    G4double Spin;

    // Coulomb Barrier
    const G4VCoulombBarrier * theCoulombBarrierPtr;

    
    // Resonances Energy
    std::vector<G4double> * ExcitationEnergies;
    
    // Resonances Spin 
    std::vector<G4double> * ExcitationSpins;

    // Resonances half lifetime
    std::vector<G4double> * ExcitationLifetimes;


    // Normalization
    G4double Normalization;
    
};


#endif
