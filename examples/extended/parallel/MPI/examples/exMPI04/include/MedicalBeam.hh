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
/// @file MedicalBeam.hh
/// @brief Define beam profile as primary generator

#ifndef MEDICAL_BEAM_H
#define MEDICAL_BEAM_H

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleDefinition;

class MedicalBeam : public G4VUserPrimaryGeneratorAction {
public:
  enum FieldShape{ kSQUARE=0, kCIRCLE };

  MedicalBeam();
  ~MedicalBeam();

  // set/get functions...
  void SetParticleDefinition(G4ParticleDefinition* pd);
  const G4ParticleDefinition* GetParticleDefinition() const;

  void SetKineticE(G4double e);
  G4double GetKineticE() const;

  void SetSourcePosition(const G4ThreeVector& pos);
  G4ThreeVector GetSourcePosition() const;

  void SetFieldShape(FieldShape shape);
  FieldShape GetFieldShape() const;

  void SetSSD(G4double ssd);
  G4double GetSSD() const;

  void SetFieldXY(G4double fx, G4double fy);
  G4double GetFieldX() const;
  G4double GetFieldY() const;

  void SetFieldR(G4double r);
  G4double GetFieldR() const;

  // methods...
  virtual void GeneratePrimaries(G4Event* anEvent);

private:
  // local methods...
  G4ThreeVector GenerateBeamDirection() const;

  G4ParticleDefinition* fparticle;
  G4double fkineticE;
  G4ThreeVector fsourcePosition;

  G4double fSSD; // (SSD= Source Skin Depth)
  FieldShape ffieldShape;
  G4double ffieldXY[2];
  G4double ffieldR;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void MedicalBeam::SetParticleDefinition(G4ParticleDefinition* pd)
{
  fparticle = pd;
}

inline const G4ParticleDefinition* MedicalBeam::GetParticleDefinition() const
{
  return fparticle;
}

inline void MedicalBeam::SetKineticE(G4double e)
{
  fkineticE = e;
}

inline G4double MedicalBeam::GetKineticE() const
{
  return fkineticE;
}

inline void MedicalBeam::SetSourcePosition(const G4ThreeVector& pos)
{
  fsourcePosition = pos;
}

inline G4ThreeVector MedicalBeam::GetSourcePosition() const
{
 return fsourcePosition;
}

inline void MedicalBeam::SetFieldShape(MedicalBeam::FieldShape shape)
{
  ffieldShape = shape;
}

inline MedicalBeam::FieldShape MedicalBeam::GetFieldShape() const
{
 return ffieldShape;
}

inline void MedicalBeam::SetSSD(G4double ssd)
{
  fSSD = ssd;
}

inline G4double MedicalBeam::GetSSD() const
{
  return fSSD;
}

inline void MedicalBeam::SetFieldXY(G4double fx, G4double fy)
{
  ffieldXY[0] = fx;
  ffieldXY[1] = fy;
}

inline G4double MedicalBeam::GetFieldX() const
{
  return ffieldXY[0];
}

inline G4double MedicalBeam::GetFieldY() const
{
  return ffieldXY[1];
}

inline void MedicalBeam::SetFieldR(G4double r)
{
  ffieldR = r;
}

inline G4double MedicalBeam::GetFieldR() const
{
  return ffieldR;
}

#endif
