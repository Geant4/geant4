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
//      File name:     G4VAngularDistribution
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 2000
//
//      Modifications: 
//      
// Abstract class for angular distribution strategy pattern
//
// Id: G4VAngularDistribution.hh,v 1.16 2000/05/11 19:07:29 pia Exp $ //
//
// -------------------------------------------------------------------

#ifndef G4VANGULARDISTRIBUTION_HH
#define G4VANGULARDISTRIBUTION_HH

#include "globals.hh"

class G4VAngularDistribution 
{

public:

  // Constructors
  G4VAngularDistribution() { }

  virtual ~G4VAngularDistribution() { }

  virtual G4double CosTheta(G4double s, G4double m1, G4double m2) const = 0;

  virtual G4double Phi() const;
 

protected:

private:  

};

#endif














