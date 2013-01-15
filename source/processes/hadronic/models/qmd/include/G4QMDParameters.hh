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
//      File name:    G4QMDParameters.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 12 May 2007
// -----------------------------------------------------------------------------

#ifndef G4QMDParameters_hh
#define G4QMDParameters_hh

#include "globals.hh"

class G4QMDParameters 
{
      static G4ThreadLocal G4QMDParameters* parameters;
     
      G4QMDParameters();
   public:
      ~G4QMDParameters();
      static G4QMDParameters* GetInstance()
      {
         if ( parameters == NULL ) parameters = new G4QMDParameters(); 
         return parameters;
      }


      // GroundStateNucleus
      G4double Get_wl() { return wl; };
      G4double Get_cl() { return cl; };
      G4double Get_hbc() { return hbc; };
      G4double Get_rho0() { return rho0; };
      G4double Get_gamm() { return gamm; };
      G4double Get_cpw() { return cpw; };
      G4double Get_cph() { return cph; };
      G4double Get_epsx() { return epsx; };
      G4double Get_cpc() { return cpc; };
      G4double Get_cs() { return cs; };
      G4double Get_c0() { return c0; };
      G4double Get_c0p() { return c0p; };
      G4double Get_clp() { return clp; };
      G4double Get_c3() { return c3; };
      G4double Get_c3p() { return c3p; };
      G4double Get_csp() { return csp; };
      G4double Get_cdp() { return cdp; };

   protected:
      G4double wl;
      G4double cl;
      G4double hbc; 
      G4double rho0;
      G4double gamm; 

/*
      G4double rho0;
      G4double t0;

      G4double aaa;

*/
      // MeanField 
      G4double c0, c3, cs;

      // GroundStateNucleus
      G4double cpw; 
      G4double cph; 
      G4double epsx;
      G4double cpc;
      G4double c0p , clp , c3p , csp , cdp;

};

#endif
