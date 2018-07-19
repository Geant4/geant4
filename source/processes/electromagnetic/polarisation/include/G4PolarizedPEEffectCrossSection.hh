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
// $Id: G4PolarizedPEEffectCrossSection.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedPEEffectCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 15.03.2007
//
// Modifications:
//
//
// Class Description:
//
#ifndef G4PolarizedPEEffectCrossSection_h
#define G4PolarizedPEEffectCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"
//#include "RotationMatrix.hh"


class G4PolarizedPEEffectCrossSection : public G4VPolarizedCrossSection
{
public:
  G4PolarizedPEEffectCrossSection();
  virtual ~G4PolarizedPEEffectCrossSection();
public:
  virtual void Initialize(G4double aGammaE, G4double aLept0E, G4double sintheta,
			  const G4StokesVector & beamPol,
			  const G4StokesVector & ,
			  G4int flag=0) override;

  G4double XSection(const G4StokesVector & pol2,
                    const G4StokesVector & pol3) override;

public:
  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;
private:
  G4StokesVector theFinalElectronPolarization;
};


#endif
