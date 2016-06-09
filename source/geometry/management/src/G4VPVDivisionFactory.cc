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
// $Id: G4VPVDivisionFactory.cc,v 1.1 2004/05/13 14:53:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// class G4VPVDivisionFactory Implementation file
//
// Author: Ivana Hrivnacova, 4.5.2004  (Ivana.Hrivnacova@cern.ch)
// ------------------------------------------------------------------------

#include "G4VPVDivisionFactory.hh"

G4VPVDivisionFactory* G4VPVDivisionFactory::fgInstance = 0;

//_____________________________________________________________________________

G4VPVDivisionFactory* G4VPVDivisionFactory::Instance() 
{
  // Static singleton access method.
  // ---
  return fgInstance;
}  


//_____________________________________________________________________________

G4VPVDivisionFactory::G4VPVDivisionFactory()
{
  // Protected singleton constructor.
  // ---
  // fgInstance = this;
}

//_____________________________________________________________________________

G4VPVDivisionFactory::~G4VPVDivisionFactory()
{
}
