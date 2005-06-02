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
// $Id: G4DNAElectronElasticEmfietzoglou.cc,v 1.1 2005-06-02 15:02:54 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DNAElectronElasticEmfietzoglou.hh"
 
                                         G4DNAElectronElasticEmfietzoglou :: G4DNAElectronElasticEmfietzoglou(const G4String & name)
:
 G4VDNAElectronElasticScatteringInWater(name),
 lowEnergyLimit(200*eV),
 highEnergyLimit(10*keV)
{
}

G4double                                 G4DNAElectronElasticEmfietzoglou :: TotalCrossSection(G4double k, G4int z)
{
 if (k<=highEnergyLimit && k>lowEnergyLimit) return (pi * RutherfordTotalCrossSection(k, z)) / ( ScreeningFactor(k,z)*(ScreeningFactor(k,z)+1.) );
 else return 0;
}


G4double                                 G4DNAElectronElasticEmfietzoglou :: RandomizeCosTheta(G4double k, G4int z)
{
 //  d sigma_el                sigma_Ruth(K)
 // ------------ (K) ~ -----------------------------  
 //   d Omega           (1 + 2 n(K) - cos(theta))^2
 //
 // We extract cos(theta) distributed as (1 + 2 n(K) - cos(theta))^-2
 //
 // Maximum is for theta=0: 1/(4 n(K)^2) (When n(K) is positive, that is always satisfied within the validity of the process)
 //
 // Phys. Med. Biol. 45 (2000) 3171-3194
 
 G4double n;
 n=ScreeningFactor(k, z);

 G4double oneOverMax;
 oneOverMax=(4.*n*n);
 
 G4double cosTheta;
 G4double fCosTheta;
 
 do
 {
  cosTheta = 2.*G4UniformRand()-1.;
  fCosTheta = (1 + 2.*n - cosTheta);
  fCosTheta = oneOverMax/(fCosTheta*fCosTheta);
 }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta;
}

