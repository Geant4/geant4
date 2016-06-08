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
//
// $Id: G4KineticTrack.hh,v 1.0 1998/05/20
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 20th May 1998
// -----------------------------------------------------------------------------

#ifndef G4KineticTrack_h
#define G4KineticTrack_h 1

#include "globals.hh"
#include "G4ios.hh"


#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4VKineticNucleon.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDecayChannel.hh"

class G4KineticTrackVector;





class G4KineticTrack : public G4VKineticNucleon
{
  public:
      
      G4KineticTrack();

      G4KineticTrack(const G4KineticTrack& right);

      G4KineticTrack(G4ParticleDefinition* aDefinition,
                     G4double aFormationTime,
                     G4ThreeVector aPosition, 
                     G4LorentzVector& a4Momentum);

      ~G4KineticTrack();

      const G4KineticTrack& operator=(const G4KineticTrack& right);

      G4int operator==(const G4KineticTrack& right) const;

      G4int operator!=(const G4KineticTrack& right) const;

      G4ParticleDefinition* GetDefinition() const;
      void SetDefinition(G4ParticleDefinition* aDefinition);

      G4double GetFormationTime() const;
      void SetFormationTime(G4double aFormationTime);

      const G4ThreeVector& GetPosition() const;
      void SetPosition(const G4ThreeVector aPosition);
      
      const G4LorentzVector& Get4Momentum() const;
      void Set4Momentum(const G4LorentzVector& a4Momentum);
      void Update4Momentum(G4double aEnergy);			// update E and p, not changing mass
      void Update4Momentum(const G4ThreeVector & aMomentum);	// idem

      const G4LorentzVector& GetInitialCoordinates() const;

      G4double SampleResidualLifetime();

      G4KineticTrackVector* Decay();
     
  // LB move to public (before was private) LB
      G4double* GetActualWidth() const;

      G4double GetActualMass() const;
      G4int GetnChannels() const;
      
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
      
  // LB new variable created LB
      G4int chosench;


  private:
 
      G4ParticleDefinition* theDefinition;

      G4double theFormationTime;

      G4ThreeVector thePosition;

      G4LorentzVector the4Momentum;
      
      G4LorentzVector theInitialCoordinates;

      G4int nChannels;
      
      G4double theActualMass;
            
      G4double* theActualWidth;



  private:

      // Temporary storage for daughter masses and widths
      // (needed because Integrand Function cannot take > 1 argument)

      G4double* theDaughterMass;

      G4double* theDaughterWidth;
};



// Class G4KineticTrack 

inline G4ParticleDefinition* G4KineticTrack::GetDefinition() const
{
  return theDefinition;
}

inline void G4KineticTrack::SetDefinition(G4ParticleDefinition* aDefinition)
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
  return the4Momentum;
}

inline void G4KineticTrack::Set4Momentum(const G4LorentzVector& a4Momentum)
{
  the4Momentum = a4Momentum;
}

inline void G4KineticTrack::Update4Momentum(G4double aEnergy)
{
  G4double newP=0;
  if ( aEnergy > GetActualMass() )
  {
      newP = sqrt(sqr(aEnergy) - the4Momentum.mag2());
  } else
  {
      aEnergy=GetActualMass();
  }
  the4Momentum = G4LorentzVector(newP*the4Momentum.vect().unit(), aEnergy);
}

inline void G4KineticTrack::Update4Momentum(const G4ThreeVector & aMomentum)
{
  G4double newE=sqrt(the4Momentum.mag2() + aMomentum.mag2());
  the4Momentum = G4LorentzVector(aMomentum, newE);
}



inline const G4LorentzVector& G4KineticTrack::GetInitialCoordinates() const
{
  return theInitialCoordinates;
}



inline G4double G4KineticTrack::GetActualMass() const
{
 G4ThreeVector theMomentum = the4Momentum.vect();
 G4double      theMomentum2 = theMomentum.mag2();
 G4double      theTotalEnergy = the4Momentum.e();
 G4double      theMass = sqrt(abs(theTotalEnergy * theTotalEnergy - theMomentum2));
 return        theMass;
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
 G4double tau = hbar_Planck * (-1.0 / theTotalActualWidth);
 G4double theResidualLifetime = tau * log(G4UniformRand());
 return theResidualLifetime*the4Momentum.gamma();
}



inline G4double G4KineticTrack::EvaluateCMMomentum(const G4double m, 
                                                 const G4double* m_ij) const
{
  G4double theCMMomentum;
  if((m_ij[0]+m_ij[1])<m)
   theCMMomentum = 1 / (2 * m) * 
          sqrt (((m * m) - (m_ij[0] + m_ij[1]) * (m_ij[0] + m_ij[1])) *
                ((m * m) - (m_ij[0] - m_ij[1]) * (m_ij[0] - m_ij[1])));
  else
   theCMMomentum=0.;

 return theCMMomentum;
}     

inline G4double G4KineticTrack::BrWig(const G4double Gamma, const G4double rmass, const G4double mass) const 
{                
  G4double Norm = twopi;
  return (Gamma/((mass-rmass)*(mass-rmass)+Gamma*Gamma/4.))/Norm;
}
#endif



