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
// $Id: G4VSDFilter.hh,v 1.1 2005/09/22 22:21:36 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

