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
// original by H.P. Wellisch
// modified by J.L. Chuma, TRIUMF, 19-Nov-1996
// last modified: 27-Mar-1997
// Chr. Volcker, 10-Nov-1997: new methods and class variables.
// M.G. Pia, 2 Oct 1998: modified GetFermiMomentum (original design was
//                       the source of memory leaks)
// G.Folger, spring 2010:  add integer A/Z interface
// A. Ribon, autumn 2021:  extended to hypernuclei

#ifndef G4Nucleus_h
#define G4Nucleus_h 1
// Class Description
// This class knows how to describe a nucleus; 
// to be used in your physics implementation (not physics list) in case you need this physics.
// Class Description - End

 
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionProductVector.hh"
#include "Randomize.hh"
 
class G4Nucleus 
{
  public:
    
    G4Nucleus();
    G4Nucleus(const G4double A, const G4double Z, const G4int numberOfLambdas = 0);
    G4Nucleus(const G4int A, const G4int Z, const G4int numberOfLambdas = 0);
    G4Nucleus(const G4Material* aMaterial);
    
    ~G4Nucleus();
    
    inline G4Nucleus( const G4Nucleus &right )
    { *this = right; }
    
    inline G4Nucleus& operator = (const G4Nucleus& right)
    {
      if (this != &right) {
        theA=right.theA;
        theZ=right.theZ;
        theL=right.theL;
        aEff=right.aEff;
        zEff=right.zEff;
        fIsotope = right.fIsotope;
        pnBlackTrackEnergy=right.pnBlackTrackEnergy; 
        dtaBlackTrackEnergy=right.dtaBlackTrackEnergy;
        pnBlackTrackEnergyfromAnnihilation =
                     right.pnBlackTrackEnergyfromAnnihilation; 
        dtaBlackTrackEnergyfromAnnihilation =
                     right.dtaBlackTrackEnergyfromAnnihilation; 
        theTemp = right.theTemp;
        excitationEnergy = right.excitationEnergy;
        momentum = right.momentum;
        fermiMomentum = right.fermiMomentum;
      }
      return *this;
    }
   
    inline G4bool operator==( const G4Nucleus &right ) const
    { return ( this == (G4Nucleus *) &right ); }
    
    inline G4bool operator!=( const G4Nucleus &right ) const
    { return ( this != (G4Nucleus *) &right ); }
    
    void ChooseParameters( const G4Material *aMaterial );

    void SetParameters( const G4double A, const G4double Z, const G4int numberOfLambdas = 0 );
    void SetParameters( const G4int A, const G4int Z, const G4int numberOfLambdas = 0 );
   
    inline G4int GetA_asInt() const
    { return theA; }   
    
    inline G4int GetN_asInt() const
    { return theA-theZ-theL; }   
    
    inline G4int GetZ_asInt() const
    { return theZ; }   

    inline G4int GetL() const  // Number of Lambdas (in the case of a hypernucleus)
    { return theL; }

    inline const G4Isotope* GetIsotope()
    { return fIsotope; }

    inline void SetIsotope(const G4Isotope* iso)
    { 
      fIsotope = iso;
      if(iso) { 
	theZ = iso->GetZ();
        theA = iso->GetN();
	theL = 0;
        aEff = theA;
        zEff = theZ;
      }
    }

    G4DynamicParticle *ReturnTargetParticle() const;
    
    G4double AtomicMass( const G4double A, const G4double Z, const G4int numberOfLambdas = 0 ) const;
    G4double AtomicMass( const G4int A, const G4int Z, const G4int numberOfLambdas = 0 ) const;

    G4double GetThermalPz( const G4double mass, const G4double temp ) const;
    
    G4ReactionProduct GetThermalNucleus(G4double aMass, G4double temp=-1) const;
    
    G4ReactionProduct GetBiasedThermalNucleus(G4double aMass, G4ThreeVector aVelocity, G4double temp=-1) const;

    G4double Cinema( G4double kineticEnergy );
    
    G4double EvaporationEffects( G4double kineticEnergy );

    G4double AnnihilationEvaporationEffects(G4double kineticEnergy, G4double ekOrg);
    
    inline G4double GetPNBlackTrackEnergy() const
    { return pnBlackTrackEnergy; }
    
    inline G4double GetDTABlackTrackEnergy() const
    { return dtaBlackTrackEnergy; }
    
    inline G4double GetAnnihilationPNBlackTrackEnergy() const
    { return pnBlackTrackEnergyfromAnnihilation; }
    
    inline G4double GetAnnihilationDTABlackTrackEnergy() const
    { return dtaBlackTrackEnergyfromAnnihilation; }
    
// ******************  methods introduced by ChV ***********************    
   // return fermi momentum
     G4ThreeVector GetFermiMomentum();

/*
  // return particle to be absorbed. 
     G4DynamicParticle* ReturnAbsorbingParticle(G4double weight);
*/

  //  final nucleus fragmentation. Return List of particles
  // which should be used for further tracking.
     G4ReactionProductVector* Fragmentate();
     

  // excitation Energy...
     void AddExcitationEnergy(G4double anEnergy);
  
  
  // momentum of absorbed Particles ..
     void AddMomentum(const G4ThreeVector aMomentum);
     
  // return excitation Energy
     G4double GetEnergyDeposit() {return excitationEnergy; }
     


// ****************************** end ChV ******************************


 private:
    
    G4int    theA;
    G4int    theZ;
    G4int    theL;  // Number of Lambdas (in the case of hypernucleus)
    G4double aEff;  // effective atomic weight
    G4double zEff;  // effective atomic number

    const G4Isotope* fIsotope;
    
    G4double pnBlackTrackEnergy;  // the kinetic energy available for
                                  // proton/neutron black track particles
    G4double dtaBlackTrackEnergy; // the kinetic energy available for
                                  // deuteron/triton/alpha particles
    G4double pnBlackTrackEnergyfromAnnihilation;
                     // kinetic energy available for proton/neutron black 
                     // track particles based on baryon annihilation 
    G4double dtaBlackTrackEnergyfromAnnihilation;
                     // kinetic energy available for deuteron/triton/alpha 
                     // black track particles based on baryon annihilation 


// ************************** member variables by ChV *******************
  // Excitation Energy leading to evaporation or deexcitation.
     G4double  excitationEnergy;
     
  // Momentum, accumulated by absorbing Particles
     G4ThreeVector momentum;
     
  // Fermi Gas model: at present, we assume constant nucleon density for all 
  // nuclei. The radius of a nucleon is taken to be 1 fm.
  // see for example S.Fl"ugge, Encyclopedia of Physics, Vol XXXIX, 
  // Structure of Atomic Nuclei (Berlin-Gottingen-Heidelberg, 1957) page 426.

  // maximum momentum possible from fermi gas model:
     G4double fermiMomentum; 
     G4double theTemp; // temperature
// ****************************** end ChV ******************************

 };
 
#endif
 
