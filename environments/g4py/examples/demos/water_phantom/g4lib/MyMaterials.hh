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
// $Id: MyMaterials.hh,v 1.2 2006-06-04 21:37:25 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyMaterials.hh
//
//                                         2005 Q
// ====================================================================
#ifndef MY_MATERIALS_H
#define MY_MATERIALS_H

#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class MyMaterials {
public:
  MyMaterials();
  ~MyMaterials();

  void Construct();
};

#endif

