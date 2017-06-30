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
// $Id: G4WentzelVIRelXSection.hh 104307 2017-05-24 09:01:45Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4WentzelVIRelXSection
//
// Authors:       V.Ivanchenko  
//
// Creation date: 08.06.2012 from G4WentzelOKandVIxSection 
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the computation of total and transport cross sections,
// sample scattering angle for the single scattering case.
// to be used by single and multiple scattering models. References:
// 1) G.Wentzel, Z. Phys. 40 (1927) 590.
// 2) J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
//
// -------------------------------------------------------------------
//

#ifndef G4WentzelVIRelXSection_h
#define G4WentzelVIRelXSection_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelOKandVIxSection.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4WentzelVIRelXSection : public G4WentzelOKandVIxSection
{

public:

  explicit G4WentzelVIRelXSection();

  virtual ~G4WentzelVIRelXSection();

  // return cos(ThetaMax) for msc and cos(thetaMin) for single scattering
  // cut = DBL_MAX means no scattering off electrons 
  virtual G4double SetupKinematic(G4double kinEnergy, const G4Material* mat);

private:

  //  hide assignment operator
  G4WentzelVIRelXSection & operator=
  (const G4WentzelVIRelXSection &right) = delete;
  G4WentzelVIRelXSection(const  G4WentzelVIRelXSection&) = delete;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

