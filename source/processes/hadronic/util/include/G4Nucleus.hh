// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Nucleus.hh,v 1.3 2000-12-14 08:56:46 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 19-Nov-1996
 // last modified: 27-Mar-1997
 // Chr. Volcker, 10-Nov-1997: new methods and class variables.
// M.G. Pia, 2 Oct 1998: modified GetFermiMomentum (original design was
//                       the source of memory leaks)
 
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
    
    G4Nucleus() { pnBlackTrackEnergy = dtaBlackTrackEnergy = 0.0;
                    excitationEnergy = 0.0;
                    momentum = G4ThreeVector(0.,0.,0.);
                    fermiMomentum = 1.52*hbarc/fermi;
                    theTemp = 293.16*kelvin;
                  }
    
    G4Nucleus( const G4double A, const G4double Z )
    {
      SetParameters( A, Z );
      pnBlackTrackEnergy = dtaBlackTrackEnergy = 0.0;
      excitationEnergy = 0.0;
      momentum = G4ThreeVector(0.,0.,0.);
      fermiMomentum = 1.52*hbarc/fermi;
      theTemp = 293.16*kelvin;
    }

    G4Nucleus( const G4Material *aMaterial )
    {
      ChooseParameters( aMaterial );
      pnBlackTrackEnergy = dtaBlackTrackEnergy = 0.0;
      excitationEnergy = 0.0;
      momentum = G4ThreeVector(0.,0.,0.);
      fermiMomentum = 1.52*hbarc/fermi;
      theTemp = aMaterial->GetTemperature();
    }
    
    ~G4Nucleus() {}
    
    inline G4Nucleus( const G4Nucleus &right )
    { *this = right; }
    
    inline G4Nucleus & operator=( const G4Nucleus &right )
     {
       if( this != &right )
       {
         aEff=right.aEff;  
         zEff=right.zEff;  
         pnBlackTrackEnergy=right.pnBlackTrackEnergy; 
         dtaBlackTrackEnergy=right.dtaBlackTrackEnergy; 
         theTemp = right.theTemp;
       }
       return *this;
     }
    
    inline G4bool operator==( const G4Nucleus &right ) const
    { return ( this == (G4Nucleus *) &right ); }
    
    inline G4bool operator!=( const G4Nucleus &right ) const
    { return ( this != (G4Nucleus *) &right ); }
    
    void ChooseParameters( const G4Material *aMaterial );

    void SetParameters( const G4double A, const G4double Z );
    
    inline G4double GetN() const
    { return aEff; }
    
    inline G4double GetZ() const
    { return zEff; }
    
    G4DynamicParticle *ReturnTargetParticle() const;
    
    G4double AtomicMass( const G4double A, const G4double Z ) const;
    
    G4double GetThermalPz( const G4double mass, const G4double temp ) const;
    
    G4ReactionProduct GetThermalNucleus(G4double aMass) const;

    G4double Cinema( G4double kineticEnergy );
    
    G4double EvaporationEffects( G4double kineticEnergy );
    
    inline G4double GetPNBlackTrackEnergy() const
    { return pnBlackTrackEnergy; }
    
    inline G4double GetDTABlackTrackEnergy() const
    { return dtaBlackTrackEnergy; }
    
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
    
    G4double aEff;  // effective atomic weight
    G4double zEff;  // effective atomic number
    
    G4double pnBlackTrackEnergy;  // the kinetic energy available for
    //                               proton/neutron black track particles
    G4double dtaBlackTrackEnergy; // the kinetic energy available for
    //                               deuteron/triton/alpha particles


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
 
