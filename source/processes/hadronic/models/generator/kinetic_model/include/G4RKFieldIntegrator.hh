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
#ifndef G4RKFieldIntegrator_h
#define G4RKFieldIntegrator_h 1

#include "G4FieldPropagation.hh"

class G4RKFieldIntegrator : public G4FieldPropagation  
{
public:
   G4RKFieldIntegrator() {}
   G4RKFieldIntegrator(const G4RKFieldIntegrator &right) {}

   ~G4RKFieldIntegrator() {}

   //Operators
   const G4RKFieldIntegrator & operator=(const G4RKFieldIntegrator &right) {return *this;}
   
   int operator==(const G4RKFieldIntegrator &right) const {return 1;}
   int operator!=(const G4RKFieldIntegrator &right) const {return 1;}

   // only theActive are propagated, nothing else
   // only theSpectators define the field, nothing else
   void Transport(G4KineticTrackVector &theActive, const G4KineticTrackVector &theSpectators, G4double theTimeStep);
   G4double GetExcitationEnergy(G4int nHitNucleons, const G4KineticTrackVector &theParticles);

   // methods for calculating potentials for different types of particles
   void Init(G4int z, G4int a) {theZ = z;  theA = a;} // prepare potentials' functions
   
   // aPosition is relative to the nucleus center
   G4double GetNeutronPotential(G4double radius);
   G4double GetNeutronPotential(G4ThreeVector &aPosition) {return GetNeutronPotential(aPosition.mag());}
   
   G4double GetProtonPotential(G4double radius);
   G4double GetProtonPotential(G4ThreeVector &aPosition)  {return GetProtonPotential(aPosition.mag());}
   
   G4double GetAntiprotonPotential(G4double radius);
   G4double GetAntiprotonPotential(G4ThreeVector &aPosition) {return GetAntiprotonPotential(aPosition.mag());};

   G4double GetKaonPotential(G4double radius);
   G4double GetKaonPotential(G4ThreeVector &aPosition) {return GetKaonPotential(aPosition.mag());}
   
   G4double GetPionPotential(G4double radius);
   G4double GetPionPotential(G4ThreeVector &aPosition) {return GetPionPotential(aPosition.mag());}
      
private:
   void Integrate(const G4KineticTrackVector & theActive, G4double theTimeStep);
   G4double CalculateTotalEnergy(const G4KineticTrackVector& Barions);
   G4double Erf(G4double X);

  // parameters to calculate potentials
  G4int theA;
  G4int theZ;
  
  // Vc(A, Z) = 1.44 * Z /(r0*(1 + pow(A, 1/3))) 
  //          = colomb * Z / (1 + pow(A, 1/3))
  static const G4double coulomb;       // coulomb barier constant
  static const G4double a_kaon;        // kaon's potential constant
  static const G4double a_pion;        // pion's potential constant
  static const G4double a_antiproton;  // antiproton's potential constant
};

#endif // G4RKFieldIntegrator_h


