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
// $Id: MaterialPropertiesGroup.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef MATERIALPROPERTIESGROUP_H
#define MATERIALPROPERTIESGROUP_H 1

#include "TagorTagref.hh"


// Morphism pattern used:
// element group -> struct
// choice        -> ContentChoice      (TagorTagref)
// choice?       -> ContentChoice* ptr (TagorTagref* ptr)
struct  MaterialPropertiesGroup
{
  // optional RL or RLref
  TagorTagref* m_RLorRLref;
  // optional AL or ALref
  TagorTagref* m_ALorALref;
  // optional T or Tref
  TagorTagref* m_TorTref;
  // optional P or Pref
  TagorTagref* m_PorPref;

  MaterialPropertiesGroup()
  : m_RLorRLref(0), m_ALorALref(0), m_TorTref(0), m_PorPref(0) {
  }
};

#endif // MATERIALPROPERTIESGROUP_H
