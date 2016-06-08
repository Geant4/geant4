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
#include "g4std/vector"

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

  G4std::vector<G4double> GenerateIso3(G4int isoIn1, G4int iso3In1, 
				       G4int isoIn2, G4int iso3In2, 
				       G4int isoOut1,G4int isoOut2) const;

  G4double Weight(G4int isoIn1,  G4int iso3In1, 
		  G4int isoIn2,  G4int iso3In2, 
		  G4int isoOut1, G4int isoOut2) const;
  
  G4double Wigner3J(G4double j1, G4double j2, G4double j3, 
		    G4double m1, G4double m2, G4double m3) const;
  
  const G4std::vector<G4double>& GetLogs() const;
  
  // Calculates the normalized Clebsch-Gordan coefficient, that is the prob 
  // of isospin decomposition of (J,m) into J1, J2, m1, m2
  G4double NormalizedClebschGordan(G4int J, G4int m, 
				   G4int J1, G4int J2,
				   G4int m1, G4int m2) const;
  
protected:

private:  

  G4Clebsch(const G4Clebsch &right);
  G4Clebsch& operator=(const G4Clebsch &right);  

  G4std::vector<G4double> logs;

};

#endif


















