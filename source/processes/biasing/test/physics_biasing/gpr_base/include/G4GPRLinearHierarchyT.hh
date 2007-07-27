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
// $Id: G4GPRLinearHierarchyT.hh,v 1.1 2007-07-27 22:13:08 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007
//
#ifndef G4GPRLINEARHIERARCHYT_HH
#define G4GPRLINEARHIERARCHYT_HH

#include "G4GPRTypeList.hh"

template <typename TList> struct G4GPRLinearHierarchyT {};

template <typename A1>
struct G4GPRLinearHierarchyT<G4GPRTypeList_1(A1)> : public A1 
{
  enum {Size = 1};
};

template <typename A1, typename A2>
struct G4GPRLinearHierarchyT<G4GPRTypeList_2(A1, A2)> : public A1, public A2 
{
  enum {Size = 2};
};

#endif
