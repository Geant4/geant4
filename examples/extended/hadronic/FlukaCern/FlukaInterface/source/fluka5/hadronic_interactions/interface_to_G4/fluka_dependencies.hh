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
// Gathers all dependencies to FLUKA,
// as needed by the G4 <-> FLUKA interface.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_DEPENDENCIES_HH
#define FLUKA_DEPENDENCIES_HH


// WRAPPERS TO FLUKA INCLUDES
#include "fluka_common_dependencies.hh"

// C++ WRAPPERS TO FLUKA PROCEDURES
#include "oauxf.h"
#include "flabrt.h"

// FLUKA PROCEDURES
extern "C" {	
	extern void cmsppr_();
	extern void zeroin_();

	extern void phycrd_(G4double* what, const char* sdum);

	extern void dcdion_(G4int& ionid);
	extern void setion_(G4int& ionid);

	extern G4double exmsaz_(G4double& aa, G4double& z, const logical& lncmss, G4int& izz0);
	extern G4double amnama_(G4double& amn, const G4int& ia, const G4int& iz);
	extern G4double enknow_(G4double& a, G4double& z, G4int& izz0);

	extern G4double xeldis_(G4int& kproj, G4double& pproj, G4int& mmat, G4double& sgmdis);

	extern void evvini_(G4double* what, const char* sdum); 		

	extern void wstoap_(G4int& niso, G4int* iaiso, G4int* iziso, G4int* iiiso, 
			    logical& lrmsch, logical& lrdhlp, logical& ltrasp, 
			    const G4int& lout);

	extern void ncdcyi_(G4int& ij, G4int& ijij, G4int& imat);
	extern void ncdcyr_();

	extern void rnread_(G4int& inseed, G4int& ijklin, logical& lseedi);
	extern void rnwrit_(const G4int& ioseed);

	extern void phncvr_(G4int& kprj, G4int& mmat, G4double& ekin, G4double& pla, G4double& dedxph);
	extern void cksigi_(const G4int& mmat);
	extern void siginm_(G4int& it, G4int& mmat, G4double& ekin, G4double& poo, G4double& sine, G4double& ainl);
	extern void pphcho_(G4double& ekin, G4int& mmat, G4double& pphnsg, const G4double& lpphch);

	extern void phnvev_(G4int& kproj, G4double& pla, G4double& ekin, 
			    G4double& txx, G4double& tyy, G4double& tzz, 
			    const G4double& txxpol, const G4double& tyypol, const G4double& tzzpol, 
			    G4double& wee, 
			    G4int& mmat, G4double& sigphn, const G4double& biaine);

	extern void eventv_(G4int& ijij, G4double& pooo, G4double& ekee, 
			    G4double& txx, G4double& tyy, G4double& tzz, 
			    G4double& we, 
			    G4int& imat);

	extern void scrprx_(G4double& eke, G4double& porenu, G4double& am, G4int& iptflg, G4double& gamcms,
			    G4double& etacms, G4double& omcsmn, G4double& omcsmx, G4double& sgscor, const G4double& tzz,
			    const logical& lpzcms);
}


#endif
#endif // G4_USE_FLUKA
