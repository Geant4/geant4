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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name: G4LightIonQMDNucleus.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 3 April 2007
// -------------------------------------------------------------------
// 230309 Skyrme-QMD parameters added by Y-H. Sato and A. Haga
// 230309 Total energy evaluated by Lorentz covariant version by Y-H. Sato and A. Haga
#ifndef G4LightIonQMDNucleus_hh
#define G4LightIonQMDNucleus_hh

#include "G4QMDSystem.hh"
#include "G4LightIonQMDParameters.hh"

class G4LightIonQMDNucleus : public G4QMDSystem
{
   public:
      G4LightIonQMDNucleus();
      //virtual ~G4LightIonQMDNucleus();

      G4LorentzVector Get4Momentum();

      // Number of Nucleons (Proton or Neutron)
      G4int GetMassNumber();

      // Number of Protons
      G4int GetAtomicNumber();

      void CalEnergyAndAngularMomentumInCM();

      // rest mass from G4NucleiPropertiesTable
      G4double GetNuclearMass();

      void SetTotalPotential( G4double x ){ potentialEnergy = x; };
      G4double GetExcitationEnergy(){ return excitationEnergy; };

      G4int GetAngularMomentum(){ return jj; };

   private:

      G4double hbc;

      std::vector < G4ThreeVector > rcm, pcm;
      std::vector < G4double > es;
      G4int jj;

      G4double potentialEnergy;
      G4double excitationEnergy;
      //G4double bindingEnergy;

      //G4double kineticEnergyPerNucleon;
      //G4double bindingEnergyPerNucleon;
      //G4double potentialEnergyPerNucleon;
    
      G4double GetNuclPotential( G4int );
    
      //G4double rclds;
    
      G4double rho0;
      G4double epsx , epscl;
    
      G4double cpc;
    
      //G4int icoul, irelcr;
      G4int irelcr;
      G4double gamm, c0, c3, cs, cl, wl;
      //G4double c0w, c3w, clw, c0sw;
      G4double c0w, clw, c0sw;
    
      G4double c0g,c3g,csg,pag;
    
      G4double cpw,cph;
    
      // Skyrme-QMD
      G4double gtau0;
      G4double g0;
      G4double g0iso;
      G4double eta;
      G4double kappas;
    
      G4double pag_tau;
      G4double cg0;
      G4double cgtau0;
};

#endif
