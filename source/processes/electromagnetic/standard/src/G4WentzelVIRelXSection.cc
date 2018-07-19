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
// $Id: G4WentzelVIRelXSection.cc 104307 2017-05-24 09:01:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelVIRelXSection
//
// Author:      V.Ivanchenko 
//
// Creation date: 08.06.2012 from G4WentzelOKandVIxSection
//
// Modifications:
//
//

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIRelXSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIRelXSection::G4WentzelVIRelXSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIRelXSection::~G4WentzelVIRelXSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4WentzelVIRelXSection::SetupKinematic(G4double kinEnergy,
	                                        const G4Material* mat)
{
  if(kinEnergy != tkin || mat != currentMaterial) { 

    currentMaterial = mat;
    tkin  = kinEnergy;
    G4double momLab2  = tkin*(tkin + 2.0*mass);
	
    G4double etot = tkin + mass;
    G4double ptot = std::sqrt(momLab2);
    G4double m12  = mass*mass;

    // relativistic reduced mass from publucation
    // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
        
    //incident particle & target nucleus
    G4double Ecm = std::sqrt(m12 + targetMass*targetMass + 2.0*etot*targetMass);
    G4double mu_rel = mass*targetMass/Ecm;
    G4double momCM  = ptot*targetMass/Ecm;
    // relative system
    mom2 = momCM*momCM;
    invbeta2 = 1.0 +  mu_rel*mu_rel/mom2;

    factB = spin/invbeta2;
    factD = std::sqrt(mom2)/targetMass;
    cosTetMaxNuc = isCombined ? 
      std::max(cosThetaMax, 1.-factorA2*mat->GetIonisation()->GetInvA23()/mom2)
      : cosThetaMax;
  } 
  return cosTetMaxNuc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

