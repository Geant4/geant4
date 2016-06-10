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
// $Id: G4tgbIsotope.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbIsotope
//
// Class description:
//
// Transient class of an isotope; builds a G4Isotope.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbIsotope_h
#define G4tgbIsotope_h

#include "globals.hh"

#include <vector>
#include <string>

#include "G4tgrIsotope.hh"
#include "G4Isotope.hh"

//-------------------------------------------------------------------------- 

class G4tgbIsotope
{

 public:   // with description 

   G4tgbIsotope();
  ~G4tgbIsotope();

   G4tgbIsotope( G4tgrIsotope* hg );
     // Construct the G4tgbIsotope from the corresponding HgIsotope

   G4Isotope* BuildG4Isotope();
     // Build a G4Isotope out of this HIsotope

   const G4String& GetName() const  { return theTgrIsot->GetName(); }

 private:

   G4tgrIsotope* theTgrIsot;
   G4Isotope* theG4Isot;
};

#endif
