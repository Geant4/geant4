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
// $Id: OlapNotify.hh,v 1.2 2003/06/16 16:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// 
// --------------------------------------------------------------
// OlapNotify
//
// Register a subclass of OlapNotify with OlapManager to be
// notified when the geometry changes or when overlap are
// reported.
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapNotify_h
#define OlapNotify_h

#include <vector>

#include "OlapEventAction.hh"

class G4LogicalVolume;

class OlapNotify
{
public:
  OlapNotify();
  virtual ~OlapNotify();
  
  virtual void worldChanged(G4LogicalVolume* newWorld) = 0;
  virtual void overlaps(const std::vector<OlapInfo*> &)=0;
				  
};								  		  		  
#endif
