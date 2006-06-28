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
// File name:     RadmonHitsManager.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHitsManager.cc,v 1.2 2006-06-28 13:57:11 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonHitsManager.hh"
#include "RadmonHit.hh"

                                                RadmonHitsManager :: RadmonHitsManager()
:
 idCounter(0),
 allocator(new Allocator)
{
}


 
 
 
G4int                                           RadmonHitsManager :: ReserveIdByLabel(const G4String & label)
{
 const IdLabels::const_iterator i(idLabels.find(label));
 
 if (i!=idLabels.end())
  return i->second;
  
 G4int id(idCounter);
 idCounter++;
 
 idLabels[label]=id;
 return id;
}

RadmonHitsManager *                             RadmonHitsManager :: instance(0);
