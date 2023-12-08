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
//      File name:    G4LightIonQMDMeanField.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 29 March 2007
// -----------------------------------------------------------------------------
// 081120 Add Update
//
// 230307 Skyrme-QMD parameters added by Y-H. Sato and A. Haga
// 230307 "CalDensityProfile" and "CalChargeDensityProfile" functions added by Y-H. Sato and A. Haga
// 230307 "GetSingleEnergy" and "GetTotalEnergy" functions added by Y-H. Sato and A. Haga

#ifndef G4LightIonQMDMeanField_hh
#define G4LightIonQMDMeanField_hh

#include "G4QMDSystem.hh"
#include "G4LightIonQMDNucleus.hh"

class G4LightIonQMDMeanField 
{
   public:

      G4LightIonQMDMeanField();
      ~G4LightIonQMDMeanField() = default;

      void SetSystem ( G4QMDSystem* aSystem );
      void SetNucleus ( G4LightIonQMDNucleus* aSystem );
      inline G4QMDSystem* GetSystem () { return system; }

      void Cal2BodyQuantities();
      void Cal2BodyQuantities( G4int );

      void CalGraduate();

      G4bool IsPauliBlocked( G4int );

      G4double GetTotalPotential();
      G4double GetPotential( G4int );

      void DoPropagation( G4double );

      std::vector< G4LightIonQMDNucleus* > DoClusterJudgment();

      inline G4double GetRR2( G4int i , G4int j ) { return rr2[i][j]; }
      inline G4double GetRHA( G4int i , G4int j ) { return rha[i][j]; }
      inline G4double GetRHE( G4int i , G4int j ) { return rhe[i][j]; }
      inline G4ThreeVector GetFFr( G4int i ) { return ffr[i]; }
      inline G4ThreeVector GetFFp( G4int i ) { return ffp[i]; }

      std::vector< G4double > GetLocalDensity();
      std::vector< G4double > GetDepthOfPotential();

      void Update();

      G4double GetTotalEnergy(); // added by Y-H.S and A.H.
      G4double GetSingleEnergy( G4int ); // added by Y-H.S and A.H.

   private:

      G4double calPauliBlockingFactor( G4int );

   private:

      G4QMDSystem* system;

      G4double rclds { 3.99 };    // distance for cluster judgement

      G4double hbc , rho0;
      G4double epsx { -20.0 };    // gauss term
      G4double epscl { 0.0001 };  // coulomb term

      G4double cpc;

      //G4int icoul, irelcr;
      G4int irelcr { 1 };
      G4double gamm, c0, c3, cs, cl, wl;
      G4double g0, g0iso, gtau0, eta, kappas; //Skyrme-QMD
      //G4double c0w, c3w, clw, c0sw;
      G4double c0w, clw, c0sw;

      G4double c0g,c3g,csg,pag; 
      G4double cg0,cgtau0,pag_tau; //Skyrme-QMD

      G4double cpw,cph;
       
      // 2 Body Quantities 
      std::vector < std::vector < G4double > > rr2;    
      std::vector < std::vector < G4double > > pp2;    
      std::vector < std::vector < G4double > > rbij;    

      // Gauss 
      std::vector < std::vector < G4double > > rha;    

      // Coulomb
      std::vector < std::vector < G4double > > rhe;    
      std::vector < std::vector < G4double > > rhc;    
                                         
      std::vector < G4ThreeVector > ffr;    
      std::vector < G4ThreeVector > ffp;    
      std::vector < G4double > rh3d;    
      std::vector < G4double > rh3d_tau; // Skyrme-QMD

};

#endif
