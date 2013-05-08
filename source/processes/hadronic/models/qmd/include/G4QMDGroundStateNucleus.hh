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
//      File name: G4QMDGroundStateNucleus.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 11 May 2007
// -------------------------------------------------------------------

#ifndef G4QMDGroundStateNucleus_hh
#define G4QMDGroundStateNucleus_hh

#include "G4QMDNucleus.hh"
#include "G4QMDMeanField.hh"

#include "G4QMDParameters.hh"

class G4QMDGroundStateNucleus : public G4QMDNucleus
{
   public:
      G4QMDGroundStateNucleus();
      G4QMDGroundStateNucleus( G4int z , G4int a );
      ~G4QMDGroundStateNucleus()
      {
         rho_l.clear();  
         d_pot.clear();
      };
                             


   private:

      void packNucleons();

      void killCMMotionAndAngularM();

      G4int maxTrial;
      G4bool samplingPosition( G4int );

      G4bool samplingMomentum( G4int );

      G4double r00 , r01 , saa , rada , radb;

      G4double dsam, ddif, dsam2, ddif2; 

      G4double cdp;
      G4double c0p,clp,c3p,csp;

      G4double hbc;
      G4double gamm;
      G4double cpw;
      G4double cph;
      G4double epsx;
      G4double cpc;

      G4double rmax,rt00,radm; 

      std::vector< G4double > phase_g;

      std::vector< G4double > rho_l;
      std::vector< G4double > d_pot;
      G4double ebini;
      G4double edepth;
      G4double epse;

      
      G4QMDMeanField* meanfield;
};

#endif
