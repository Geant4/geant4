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
// $Id: G4AblaFissionSimfis18.hh,v 1.2 2008-06-25 17:20:04 kaitanie Exp $
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#ifndef G4AblaFissionSimfis18_hh
#define G4AblaFissionSimfis18_hh 1

#include "G4AblaFissionBase.hh"
#include "G4InclDataDefs.hh"
#include "G4InclRandomNumbers.hh"

class G4AblaFissionSimfis18 : public G4AblaFissionBase {

public:
  G4AblaFissionSimfis18();
  G4AblaFissionSimfis18(G4Hazard *hzr, G4InclRandomInterface *rndm);

  ~G4AblaFissionSimfis18();

  void doFission(G4double &A, G4double &Z, G4double &E,
		 G4double &A1, G4double &Z1, G4double &E1, G4double &K1,
		 G4double &A2, G4double &Z2, G4double &E2, G4double &K2);

  /**
   *
   */
  G4double spdef(G4int a, G4int z, G4int optxfis);

  /**
   *
   */
  G4double fissility(G4int a, G4int z, G4int optxfis);

//   void evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf,
// 	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
// 	       G4double *pleva_par, G4double *pxeva_par);
//  G4double bfms67(G4double zms, G4double ams);
  //  void lpoly(G4double x, G4int n, G4double pl[]);
  //  G4double expohaz(G4int k, G4double T);
  //  G4double fd(G4double E);
  //  G4double f(G4double E);
  //  G4double fmaxhaz(G4double k, G4double T);
  void even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out);
  G4double umass(G4double z,G4double n,G4double beta);
  G4double ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d);
  void fissionDistri(G4double &a,G4double &z,G4double &e,
		     G4double &a1,G4double &z1,G4double &e1,G4double &v1,
		     G4double &a2,G4double &z2,G4double &e2,G4double &v2);
  void standardRandom(G4double *rndm, G4long *seed);
  G4double haz(G4int k);
  G4double gausshaz(int k, double xmoy, double sig);



  G4int min(G4int a, G4int b);
  G4double min(G4double a, G4double b);
  G4int max(G4int a, G4int b);
  G4double max(G4double a, G4double b);

  G4int nint(G4double number);
  G4int secnds(G4int x);
  G4int mod(G4int a, G4int b);
  G4double dmod(G4double a, G4double b);
  G4double dint(G4double a);
  G4int idint(G4double a);
  G4int idnint(G4double value);
  G4double utilabs(G4double a);
  G4double dmin1(G4double a, G4double b, G4double c);

private:
  G4InclRandomInterface *randomGenerator;
  G4Hazard *hazard;
};

#endif
