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
#ifndef G4RKFieldIntegrator_h
#define G4RKFieldIntegrator_h 1

#include "G4FieldPropagation.hh"

class G4RKFieldIntegrator : public G4FieldPropagation  
{
public:
   G4RKFieldIntegrator() {}
   G4RKFieldIntegrator(const G4RKFieldIntegrator &):G4FieldPropagation() {}

   ~G4RKFieldIntegrator() {}

   //Operators
   const G4RKFieldIntegrator & operator=(const G4RKFieldIntegrator &) {return *this;}
   
   int operator==(const G4RKFieldIntegrator &) const {return 1;}
   int operator!=(const G4RKFieldIntegrator &) const {return 1;}

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
  
  // Vc(A, Z) = 1.44 * Z /(r0*(1 + std::pow(A, 1/3))) 
  //          = colomb * Z / (1 + std::pow(A, 1/3))
  static const G4double coulomb;       // coulomb barier constant
  static const G4double a_kaon;        // kaon's potential constant
  static const G4double a_pion;        // pion's potential constant
  static const G4double a_antiproton;  // antiproton's potential constant
};

#endif // G4RKFieldIntegrator_h


