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
// $Id: NTSTGunGenerator.hh,v 1.2 2003-12-09 15:35:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTGunGenerator_h
#define NTSTGunGenerator_h 1

#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"

class G4Event;
class NTSTGunMessenger;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTGunGenerator : public G4VPrimaryGenerator
{
public:
  NTSTGunGenerator();
  virtual ~NTSTGunGenerator();
  virtual void GeneratePrimaryVertex(G4Event* evt);

public:
  virtual void SetParticleDefinition (G4ParticleDefinition* pDef)
    {particle_definition = pDef;}
  virtual void SetPlow        (const G4double m)      { Plow  = m; }
  virtual void SetPhigh       (const G4double m)      { Phigh = m; }
  virtual void SetCoslow      (const G4double m)      { Coslow  = m; }
  virtual void SetCoshigh     (const G4double m)      { Coshigh = m; }
  virtual void SetPhilow      (const G4double m)      { Philow  = m; }
  virtual void SetPhihigh     (const G4double m)      { Phihigh = m; }
  virtual void SetMeanVertex  (const G4ThreeVector v) { MeanVertex = v;}
  virtual void SetRmsVertex   (const G4ThreeVector v) { RmsVertex  = v;}
  virtual void SetPolarization(const G4ThreeVector v) { Polarization= v;}
  virtual void SetT0          (const G4double t)      { T0=t; }
  virtual void SetNumberOfParticles   (const G4int n) { N=n; }

public:
  virtual G4ParticleDefinition* GetParticleDefinition()
    {return particle_definition;}
  virtual const G4ThreeVector& GetMeanVertex()   const {return MeanVertex;}
  virtual const G4ThreeVector& GetRmsVertex()    const {return RmsVertex;}
  virtual const G4ThreeVector& GetPolarization() const {return Polarization;}
  virtual const G4double       GetPlow()         const {return Plow;}
  virtual const G4double       GetPhigh()        const {return Phigh;}
  virtual const G4double       GetCoslow()       const {return Coslow;}
  virtual const G4double       GetCoshigh()      const {return Coshigh;}
  virtual const G4double       GetPhilow()       const {return Philow; }
  virtual const G4double       GetPhihigh()      const {return Phihigh;}
  virtual const G4double       GetT0()           const {return T0;}
  virtual const G4int          GetNumberOfParticles() const {return N;}

protected:
  inline  G4double Gauss(){
    return RandGauss::shoot();
  }

private:
  G4ParticleDefinition* particle_definition;
  G4ThreeVector         MeanVertex;
  G4ThreeVector         RmsVertex;
  G4ThreeVector         Polarization;
  G4double              Plow;
  G4double              Phigh;
  G4double              Coslow;
  G4double              Coshigh;
  G4double              Philow;
  G4double              Phihigh;
  G4double              T0;
  G4int                 N;

  G4ThreeVector         dMeanVertex;
  G4ThreeVector         dRmsVertex;
  G4ThreeVector         dPolarization;
  G4double              dPlow;
  G4double              dPhigh;
  G4double              dCoslow;
  G4double              dCoshigh;
  G4double              dPhilow;
  G4double              dPhihigh;
  G4double              dT0;
  G4int                 dN;

  NTSTGunMessenger*   messenger;
};

#endif


