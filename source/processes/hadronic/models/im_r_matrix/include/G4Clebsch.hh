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
//
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4Clebsch
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
//      Kinetic Model, Clebsch-Gordan coefficient and related algebra
//
// -------------------------------------------------------------------

#ifndef G4CLEBSCH_HH
#define G4CLEBSCH_HH

#include "globals.hh"
#include <vector>

class G4Clebsch 
{

public:

  G4Clebsch();

  virtual ~G4Clebsch();

  G4bool operator==(const G4Clebsch &right) const;
  G4bool operator!=(const G4Clebsch &right) const;

  G4double ClebschGordan(G4int isoIn1, G4int iso3In1, 
			 G4int isoIn2, G4int iso3In2, 
			 G4int jOut) const;

  std::vector<G4double> GenerateIso3(G4int isoIn1, G4int iso3In1, 
				       G4int isoIn2, G4int iso3In2, 
				       G4int isoOut1,G4int isoOut2) const;

  G4double Weight(G4int isoIn1,  G4int iso3In1, 
		  G4int isoIn2,  G4int iso3In2, 
		  G4int isoOut1, G4int isoOut2) const;
  
  G4double Wigner3J(G4double j1, G4double j2, G4double j3, 
		    G4double m1, G4double m2, G4double m3) const;
  
  // Calculates the normalized Clebsch-Gordan coefficient, that is the prob 
  // of isospin decomposition of (J,m) into J1, J2, m1, m2
  G4double NormalizedClebschGordan(G4int J, G4int m, 
				   G4int J1, G4int J2,
				   G4int m1, G4int m2) const;
  
protected:

private:  

  G4Clebsch(const G4Clebsch &right);
  G4Clebsch& operator=(const G4Clebsch &right);  

  std::vector<G4double> logs;

};

#endif


















