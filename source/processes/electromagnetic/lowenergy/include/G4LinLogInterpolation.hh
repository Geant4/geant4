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
//
// Author:  Vladimir Ivanchenko (Vladimir.Ivantchenko@cern.ch)
//
// History:
// -----------
// 27 May 2002   VI        Created
//
// -------------------------------------------------------------------
// Class description:
// Lin-Log interpolation of a data set

// -------------------------------------------------------------------

#ifndef G4LINLOGINTERPOLATION_HH
#define G4LINLOGINTERPOLATION_HH 1

#include "globals.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4DataVector.hh"

class G4LinLogInterpolation : public G4VDataSetAlgorithm {
 
public:

  G4LinLogInterpolation();

  ~G4LinLogInterpolation();
 
  G4double Calculate(G4double point, G4int bin, 
		     const G4DataVector& energies, 
		     const G4DataVector& data) const;

  virtual G4VDataSetAlgorithm* Clone() const { return new G4LinLogInterpolation; }

private:

  
  // Hide copy constructor and assignment operator

};
 
#endif
 










