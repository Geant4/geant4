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
// $Id: TiaraIsotropicDirections.hh,v 1.4 2006/06/29 15:43:54 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraIsotropicDirections
//

#ifndef TiaraIsotropicDirections_hh
#define TiaraIsotropicDirections_hh TiaraIsotropicDirections_hh

#include "TiaraVDirectionGenerator.hh"
#include "G4ThreeVector.hh"

class TiaraDimensions;

class TiaraIsotropicDirections : public TiaraVDirectionGenerator {
public:
  TiaraIsotropicDirections(G4double colWidth,
			   const TiaraDimensions &tiaraDimensions);
  ~TiaraIsotropicDirections();
  TiaraIsotropicDirections(const TiaraIsotropicDirections& rhs);

  virtual G4ThreeVector GetDirection();
  virtual TiaraVDirectionGenerator *Clone() const;

  G4double MinimumCosine(G4double colWidth,
			 const TiaraDimensions &tiaraDimensions);
  
  TiaraIsotropicDirections& operator=(const TiaraIsotropicDirections& rhs);
private:
  G4double fMinCos;
};

#endif
