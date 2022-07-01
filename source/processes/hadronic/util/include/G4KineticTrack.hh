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
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 20th May 1998
// -----------------------------------------------------------------------------

#ifndef G4KineticTrack_h
#define G4KineticTrack_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"


#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4VKineticNucleon.hh"
#include "G4Nucleon.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDecayChannel.hh"
#include "G4Log.hh"

// #include "G4Allocator.hh"

class G4KineticTrackVector;

class G4KineticTrack : public G4VKineticNucleon
{
  public:
      
      G4KineticTrack();

      G4KineticTrack(const G4KineticTrack& right);

      G4KineticTrack(const G4ParticleDefinition* aDefinition,
                     G4double aFormationTime,
                     const G4ThreeVector& aPosition,
                     const G4LorentzVector& a4Momentum);
      G4KineticTrack(G4Nucleon * nucleon,
                     const G4ThreeVector& aPosition,
                     const G4LorentzVector& a4Momentum);

      ~G4KineticTrack();

      G4KineticTrack& operator=(const G4KineticTrack& right);

      G4bool operator==(const G4KineticTrack& right) const;

      G4bool operator!=(const G4KineticTrack& right) const;
/*
      inline void *operator new(size_t);
      inline void operator delete(void *aTrack);
*/
      const G4ParticleDefinition* GetDefinition() const;
      void SetDefinition(const G4ParticleDefinition* aDefinition);

      G4double GetFormationTime() const;
      void SetFormationTime(G4double aFormationTime);

      const G4ThreeVector& GetPosition() const;
      void SetPosition(const G4ThreeVector aPosition);
      
      const G4LorentzVector& Get4Momentum() const;
      void Set4Momentum(const G4LorentzVector& a4Momentum);
      void Update4Momentum(G4double aEnergy);			// update E and p, not changing mass
      void Update4Momentum(const G4ThreeVector & aMomentum);	// idem
      void SetTrackingMomentum(const G4LorentzVector& a4Momentum);
      void UpdateTrackingMomentum(G4double aEnergy);			// update E and p, not changing mass
      void UpdateTrackingMomentum(const G4ThreeVector & aMomentum);	// idem

      const G4LorentzVector& GetTrackingMomentum() const;
      
      G4double SampleResidualLifetime();
      
      void Hit(); 
      void SetNucleon(G4Nucleon * aN) {theNucleon = aN;}
            
      G4bool IsParticipant() const; 

      G4KineticTrackVector* Decay();
     
  // LB move to public (before was private) LB
      G4double* GetActualWidth() const;

      G4double GetActualMass() const;
      G4int GetnChannels() const;
      
//   position relativ to nucleus "state"
      enum CascadeState {undefined, outside, going_in, inside, 
                         going_out, gone_out, captured, miss_nucleus };
      
      CascadeState SetState(const CascadeState new_state);
      CascadeState GetState() const;
      void SetProjectilePotential(const G4double aPotential);
      G4double GetProjectilePotential() const;

      void SetCreatorModelID(G4int id);
      G4int GetCreatorModelID() const;

      const G4ParticleDefinition* GetParentResonanceDef() const;
      void SetParentResonanceDef(const G4ParticleDefinition* parent);      
      G4int GetParentResonanceID() const;
      void SetParentResonanceID(const G4int parentID);
   
  private:


      void SetnChannels(const G4int aChannel);

      void SetActualWidth(G4double* anActualWidth); 
      
      G4double EvaluateTotalActualWidth();

      G4double EvaluateCMMomentum (const G4double mass,
                                   const G4double* m_ij) const;                                 
      
      G4double IntegrateCMMomentum(const G4double lowerLimit) const;

      G4double IntegrateCMMomentum(const G4double lowerLimit ,const G4double polemass) const;

      G4double IntegrateCMMomentum2() const;
      
  public:
      
      G4double BrWig(const G4double Gamma, 
                     const G4double rmass, 
                     const G4double mass) const;

private:
      G4double IntegrandFunction1 (G4double xmass) const;
      G4double IntegrandFunction2 (G4double xmass) const;
      G4double IntegrandFunction3 (G4double xmass) const;
      G4double IntegrandFunction4 (G4double xmass) const;
public:
  //   friend G4double IntegrandFunction3 (G4double xmass);

  //   friend G4double IntegrandFunction4 (G4double xmass);
      
  private:
 
      const G4ParticleDefinition* theDefinition;

      G4double theFormationTime;

      G4ThreeVector thePosition;

      G4LorentzVector the4Momentum;
      G4LorentzVector theFermi3Momentum;
      G4LorentzVector theTotal4Momentum;
      
      G4Nucleon * theNucleon;
      
      G4int nChannels;
      
      G4double theActualMass;
            
      G4double* theActualWidth;

     // Temporary storage for daughter masses and widths
      // (needed because Integrand Function cannot take > 1 argument)
      G4double* theDaughterMass;
      G4double* theDaughterWidth;

      CascadeState theStateToNucleus;

      G4double theProjectilePotential;

      G4int theCreatorModel;

      const G4ParticleDefinition* theParentResonanceDef = nullptr;
      G4int theParentResonanceID;
};

// extern G4Allocator<G4KineticTrack> theKTAllocator;


// Class G4KineticTrack 
/*
inline void * G4KineticTrack::operator new(size_t)
{
  void * aT;
  aT = (void *) theKTAllocator.MallocSingle();
  return aT;
}

inline void G4KineticTrack::operator delete(void * aT)
{
  theKTAllocator.FreeSingle((G4KineticTrack *) aT);
}
*/

inline const G4ParticleDefinition* G4KineticTrack::GetDefinition() const
{
  return theDefinition;
}

inline void G4KineticTrack::SetDefinition(const G4ParticleDefinition* aDefinition)
{
  theDefinition = aDefinition;
}

inline G4double G4KineticTrack::GetFormationTime() const
{
  return theFormationTime;
}

inline void G4KineticTrack::SetFormationTime(G4double aFormationTime)
{
  theFormationTime = aFormationTime;
}

inline const G4ThreeVector& G4KineticTrack::GetPosition() const
{
  return thePosition;
}

inline void G4KineticTrack::SetPosition(const G4ThreeVector aPosition)
{
  thePosition = aPosition;
}

inline const G4LorentzVector& G4KineticTrack::Get4Momentum() const
{
  return theTotal4Momentum;
}

inline const G4LorentzVector& G4KineticTrack::GetTrackingMomentum() const
{
   return the4Momentum;
}

inline void G4KineticTrack::Set4Momentum(const G4LorentzVector& a4Momentum)
{
//  set the4Momentum and update theTotal4Momentum

  theTotal4Momentum=a4Momentum;
  the4Momentum = theTotal4Momentum;
  theFermi3Momentum=G4LorentzVector(0);
}

inline void G4KineticTrack::Update4Momentum(G4double aEnergy)
{
// update the4Momentum with aEnergy at constant mass (the4Momentum.mag()  
//   updates theTotal4Momentum as well.
  G4double newP(0);
  G4double mass2=theTotal4Momentum.mag2();
  if ( sqr(aEnergy) > mass2 )
  {
      newP = std::sqrt(sqr(aEnergy) - mass2 );
  } else
  {
      aEnergy=std::sqrt(mass2);
  }
  Set4Momentum(G4LorentzVector(newP*the4Momentum.vect().unit(), aEnergy));
}

inline void G4KineticTrack::Update4Momentum(const G4ThreeVector & aMomentum)
{
// update the4Momentum with aMomentum at constant mass (the4Momentum.mag()  
//   updates theTotal4Momentum as well.
  G4double newE=std::sqrt(theTotal4Momentum.mag2() + aMomentum.mag2());
  Set4Momentum(G4LorentzVector(aMomentum, newE));
}

inline void G4KineticTrack::SetTrackingMomentum(const G4LorentzVector& aMomentum)
{
//  set the4Momentum and update theTotal4Momentum, keep the mass of aMomentum

  the4Momentum = aMomentum;
  theTotal4Momentum=the4Momentum+theFermi3Momentum;
//     keep mass of aMomentum for the total momentum
  G4double mass2 = aMomentum.mag2();
  G4double p2=theTotal4Momentum.vect().mag2();
  theTotal4Momentum.setE(std::sqrt(mass2+p2));
}

inline void G4KineticTrack::UpdateTrackingMomentum(G4double aEnergy)
{
// update the4Momentum with aEnergy at constant mass (the4Momentum.mag()  
//   updates theTotal4Momentum as well.
  G4double newP(0);
  G4double mass2=theTotal4Momentum.mag2();
  if ( sqr(aEnergy) > mass2 )
  {
      newP = std::sqrt(sqr(aEnergy) - mass2 );
  } else
  {
      aEnergy=std::sqrt(mass2);
  }
  SetTrackingMomentum(G4LorentzVector(newP*the4Momentum.vect().unit(), aEnergy));
}

inline void G4KineticTrack::UpdateTrackingMomentum(const G4ThreeVector & aMomentum)
{
// update the4Momentum with aMomentum at constant mass (the4Momentum.mag()  
//   updates theTotal4Momentum as well.
  G4double newE=std::sqrt(theTotal4Momentum.mag2() + aMomentum.mag2());
  SetTrackingMomentum(G4LorentzVector(aMomentum, newE));
}

inline G4double G4KineticTrack::GetActualMass() const
{
  return std::sqrt(std::abs(the4Momentum.mag2()));
}

inline G4int G4KineticTrack::GetnChannels() const
{
  return nChannels;
}

inline void G4KineticTrack::SetnChannels(const G4int numberOfChannels)
{
  nChannels = numberOfChannels;
}

inline G4double* G4KineticTrack::GetActualWidth() const
{
  return theActualWidth;
}

inline void G4KineticTrack::SetActualWidth(G4double* anActualWidth)
{
  theActualWidth = anActualWidth;
}



inline G4double G4KineticTrack::EvaluateTotalActualWidth()
{
 G4int index;
 G4double theTotalActualWidth = 0.0;
 for (index = nChannels - 1; index >= 0; index--)
    {
     theTotalActualWidth += theActualWidth[index];
    }
 return theTotalActualWidth;
}

inline G4double G4KineticTrack::SampleResidualLifetime()
{
 G4double theTotalActualWidth = this->EvaluateTotalActualWidth();
 G4double tau = CLHEP::hbar_Planck * (-1.0 / theTotalActualWidth);
 G4double theResidualLifetime = tau * G4Log(G4UniformRand());
 return theResidualLifetime*the4Momentum.gamma();
}

inline G4double G4KineticTrack::EvaluateCMMomentum(const G4double mass, 
                                                 const G4double* m_ij) const
{
  G4double theCMMomentum;
  if((m_ij[0]+m_ij[1])<mass)
   theCMMomentum = 1 / (2 * mass) * 
          std::sqrt (((mass * mass) - (m_ij[0] + m_ij[1]) * (m_ij[0] + m_ij[1])) *
                ((mass * mass) - (m_ij[0] - m_ij[1]) * (m_ij[0] - m_ij[1])));
  else
   theCMMomentum=0.;

 return theCMMomentum;
}     

inline G4double G4KineticTrack::BrWig(const G4double Gamma, const G4double rmass, const G4double mass) const 
{                
  G4double Norm = CLHEP::twopi;
  return (Gamma/((mass-rmass)*(mass-rmass)+Gamma*Gamma/4.))/Norm;
}
      
inline      
void G4KineticTrack::Hit() 
{
  if(theNucleon) 
  {
    theNucleon->Hit(1);
  }
}

inline
G4bool G4KineticTrack::IsParticipant() const 
{ 
  if(!theNucleon) return true;
  return theNucleon->AreYouHit(); 
}

inline 
G4KineticTrack::CascadeState G4KineticTrack::GetState() const
{
	return theStateToNucleus;
}

inline
G4KineticTrack::CascadeState G4KineticTrack::SetState(const CascadeState new_state)
{
	CascadeState old_state=theStateToNucleus;
	theStateToNucleus=new_state;
	return old_state;
}

inline
void G4KineticTrack::SetProjectilePotential(G4double aPotential)
{
	theProjectilePotential = aPotential;
}
inline
G4double G4KineticTrack::GetProjectilePotential() const
{
	return theProjectilePotential;
}

inline
void G4KineticTrack::SetCreatorModelID(G4int id)
{
        theCreatorModel = id;
}
inline
G4int G4KineticTrack::GetCreatorModelID() const
{
       return theCreatorModel;
}

inline
const G4ParticleDefinition* G4KineticTrack::GetParentResonanceDef() const
{
        return theParentResonanceDef;
}

inline
void G4KineticTrack::SetParentResonanceDef(const G4ParticleDefinition* parentDef)
{
        theParentResonanceDef = parentDef;
}

inline
G4int G4KineticTrack::GetParentResonanceID() const
{
        return theParentResonanceID;
}

inline
void G4KineticTrack::SetParentResonanceID(const G4int parentID)
{
        theParentResonanceID = parentID;
}

#endif
