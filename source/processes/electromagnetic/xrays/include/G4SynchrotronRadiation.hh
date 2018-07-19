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
// $Id: G4SynchrotronRadiation.hh 97385 2016-06-02 09:59:53Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      History:
//      21-5-98  1 version , V. Grichine
//      28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
//      23-05-06, H. Burkhardt: Energy spectrum from function rather than table
//
//
// ------------------------------------------------------------

#ifndef G4SynchrotronRadiation_h
#define G4SynchrotronRadiation_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"
#include "G4PropagatorInField.hh"

#include "G4Track.hh"
#include "G4Step.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4VEmAngularDistribution;

class G4SynchrotronRadiation : public G4VDiscreteProcess
{
public:

  explicit G4SynchrotronRadiation(const G4String& pName = "SynRad",
		         G4ProcessType type = fElectromagnetic);

  virtual ~G4SynchrotronRadiation();

  virtual G4double GetMeanFreePath( const G4Track& track,
				    G4double previousStepSize,
				    G4ForceCondition* condition ) override;

  virtual G4VParticleChange *PostStepDoIt( const G4Track& track,
					   const G4Step& Step    ) override;

  G4double GetPhotonEnergy( const G4Track& trackData,
                               const G4Step&  stepData      );

  G4double GetRandomEnergySR( G4double, G4double, G4double );

  G4double InvSynFracInt(G4double x);
  G4double Chebyshev(G4double a,G4double b,const G4double c[],
		     G4int n, G4double x);

  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;
  virtual void BuildPhysicsTable(const G4ParticleDefinition& ) override;
  virtual void PrintInfoDefinition();

  void SetAngularGenerator(G4VEmAngularDistribution* p);

private:

  G4SynchrotronRadiation & 
    operator=(const G4SynchrotronRadiation &right) = delete;
  G4SynchrotronRadiation(const G4SynchrotronRadiation&) = delete;

  G4VEmAngularDistribution*   genAngle;

  G4ParticleDefinition*       theGamma;

  G4PropagatorInField* fFieldPropagator;
  G4bool FirstTime;
  G4bool FirstTime1;
};

//////////////////////////  INLINE METHODS  /////////////////////////////

inline G4double 
G4SynchrotronRadiation::Chebyshev(G4double a, G4double b, const G4double c[],
				  G4int n, G4double x)
{
  G4double y;
  G4double y2=2.0*(y=(2.0*x-a-b)/(b-a)); // Change of variable.
  G4double d=0.,dd=0.;
  for (G4int j=n-1;j>=1;--j) // Clenshaw's recurrence.
  { G4double sv=d;
	d=y2*d-dd+c[j];
	dd=sv;
  }
  return y*d-dd+0.5*c[0];
}

#endif  // end of G4SynchrotronRadiation.hh

