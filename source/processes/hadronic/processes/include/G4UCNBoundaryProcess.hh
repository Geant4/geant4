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
//
////////////////////////////////////////////////////////////////////////
// Ultra Cold Neutron (UCN) Boundary Process Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:         G4UCNBoundaryProcess.hh
// Description:  Discrete Process -- Boundary Process of UCN
// Version:      1.0
// Created:      2014-06-04
// Author:       Peter Gumplinger
// Adopted from: UCNMaterialBoundary by Peter Fierlinger 4.9.2004
// Updated:
// mail:         gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4UCNBOUNDARYPROCESS_HH
#define G4UCNBOUNDARYPROCESS_HH 1

/////////////
// Includes
/////////////

#include "G4VDiscreteProcess.hh"

#include "G4Neutron.hh"

#include "G4UCNMaterialPropertiesTable.hh"

class G4UCNBoundaryProcessMessenger;

// Class Description:
// Discrete Process -- Boundary Process of Ultra Cold Neutrons.
// Reflects/Absorpts UCN at boundaries.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////


enum G4UCNBoundaryProcessStatus { Undefined,
                                  NotAtBoundary,
                                  SameMaterial,
                                  StepTooSmall,
                                  NoMPT, NoMRT,
                                  NoMRCondition,
                                  Absorption, Ezero, Flip,
                                  SpecularReflection,
                                  LambertianReflection, MRDiffuseReflection,
                                  SnellTransmit, MRDiffuseTransmit
                                };

class G4UCNBoundaryProcess : public G4VDiscreteProcess
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4UCNBoundaryProcess(const G4String& processName = "UCNBoundaryProcess",
                             G4ProcessType type = fUCN);
	virtual ~G4UCNBoundaryProcess();

private:

        G4UCNBoundaryProcess(const G4UCNBoundaryProcess &right);

        //////////////
        // Operators
        //////////////

        G4UCNBoundaryProcess& operator=(const G4UCNBoundaryProcess &right);

public:

        ////////////
        // Methods
        ///////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an UCN.

        G4double GetMeanFreePath(const G4Track& aTrack,
                                 G4double ,
                                 G4ForceCondition* condition);

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);

private:

        G4UCNBoundaryProcessMessenger* fMessenger;

        G4double neV;

        G4double kCarTolerance;

        G4UCNBoundaryProcessStatus theStatus;

        const G4Material* Material1;
        const G4Material* Material2;

        // the G4UCNMaterialPropertiesTable of PreStepPoint
        G4UCNMaterialPropertiesTable* aMaterialPropertiesTable1;
        // the G4UCNMaterialPropertiesTable of PostStepPoint
        G4UCNMaterialPropertiesTable* aMaterialPropertiesTable2;

        G4bool UseMicroRoughnessReflection;
        G4bool DoMicroRoughnessReflection;

        ////////////
        // Methods
        ///////////

        G4bool High(G4double , G4double );

        G4bool Loss(G4double , G4double , G4double );

        G4bool SpinFlip(G4double );

        G4double Reflectivity(G4double , G4double );

        G4ThreeVector Reflect(G4double , G4ThreeVector , G4ThreeVector );

        G4double Transmit(G4double, G4double );

        G4ThreeVector LDiffRefl(G4ThreeVector );

public:
  
        G4ThreeVector MRreflect(G4double ,
                                G4ThreeVector , G4ThreeVector ,
                                G4double , G4double );
        G4ThreeVector MRreflectHigh(G4double , G4double , G4double ,
                                    G4ThreeVector , G4ThreeVector ,
                                    G4double , G4double , G4double& );

private:

        G4ThreeVector MRDiffRefl(G4ThreeVector , G4double , G4double ,
                                 G4ThreeVector , G4double );
        G4ThreeVector MRDiffTrans(G4ThreeVector , G4double , G4double ,
                                  G4ThreeVector , G4double );

        G4RotationMatrix GetCoordinateTransformMatrix(G4ThreeVector ,
                                                      G4ThreeVector );

        void BoundaryProcessVerbose() const;

        // Invoke SD for post step point if the photon is 'detected'
        G4bool InvokeSD(const G4Step* step);

private:

        G4int nNoMPT, nNoMRT, nNoMRCondition;
        G4int nAbsorption, nEzero, nFlip;
        G4int aSpecularReflection, bSpecularReflection;
        G4int bLambertianReflection;
        G4int aMRDiffuseReflection, bMRDiffuseReflection;
        G4int nSnellTransmit, mSnellTransmit;
        G4int aMRDiffuseTransmit;

        G4double ftheta_o, fphi_o;

public:

        void SetMicroRoughness(G4bool );
        G4bool GetMicroRoughness();

        G4UCNBoundaryProcessStatus GetStatus() const;

        void BoundaryProcessSummary() const;

        void SetMaterialPropertiesTable1(G4UCNMaterialPropertiesTable* MPT)
             {aMaterialPropertiesTable1 = MPT;}

        void SetMaterialPropertiesTable2(G4UCNMaterialPropertiesTable* MPT)
             {aMaterialPropertiesTable2 = MPT;}

        G4double GetTheta_o() {return ftheta_o;};
        G4double GetPhi_o() {return fphi_o;};

};

////////////////////
// Inline methods
////////////////////

inline G4bool
G4UCNBoundaryProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return ( &aParticleType == G4Neutron::NeutronDefinition() );
}

inline
G4UCNBoundaryProcessStatus G4UCNBoundaryProcess::GetStatus() const
{
  return theStatus;
}

inline
G4bool G4UCNBoundaryProcess::High(G4double Energy, G4double FermiPotDiff)
{
  // Returns true for Energy > Fermi Potential Difference

  return (Energy > FermiPotDiff);
}

inline
void G4UCNBoundaryProcess::SetMicroRoughness(G4bool active)
{
   UseMicroRoughnessReflection = active;
}

inline
G4bool G4UCNBoundaryProcess::GetMicroRoughness()
{
  return UseMicroRoughnessReflection;
}

#endif /* G4UCNBOUNDARYPROCESS_HH */
