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
//      File name:     G4ResonanceWidth
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4RESONANCEWIDTH_HH
#define G4RESONANCEWIDTH_HH

#include <map>
#include "globals.hh"

class G4PhysicsVector;

class G4ResonanceWidth 
{
public:

  G4ResonanceWidth() {};

  virtual ~G4ResonanceWidth() {};

  // Returned pointer is owned by the client
  virtual G4PhysicsVector* MassDependentWidth(const G4String& name) const = 0;

protected:
  
private:  

  G4ResonanceWidth(const G4ResonanceWidth& right);
  G4ResonanceWidth& operator=(const G4ResonanceWidth& right);
};
  
#endif
