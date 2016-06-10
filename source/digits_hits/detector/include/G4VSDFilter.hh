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
// $Id: G4VSDFilter.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4VSDFilter_h
#define G4VSDFilter_h 1

class G4Step;
#include "globals.hh"

// class description:
//
//  This is the abstract base class of a filter to be associated with a
// sensitive detector. 

class G4VSDFilter 
{

  public: // with description
      G4VSDFilter(G4String name);

  public:
      virtual ~G4VSDFilter();

  public: // with description
      virtual G4bool Accept(const G4Step*) const = 0;

  protected:
      G4String filterName;

  public:
      inline G4String GetName() const
      { return filterName; }
};

#endif

