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
// $Id: G4GPRLinearHierarchyT.hh,v 1.4 2007-09-14 16:42:50 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007
//
#ifndef G4GPRLINEARHIERARCHYT_HH
#define G4GPRLINEARHIERARCHYT_HH

#include "G4GPRTypeList.hh"

template <typename TList> struct G4GPRLinearHierarchyT {};

template <typename A1>
struct G4GPRLinearHierarchyT<G4GPRTypeList_1(A1)> : public A1 {};

template <typename A1, typename A2>
struct G4GPRLinearHierarchyT<G4GPRTypeList_2(A1, A2)> : public A1, public A2 {};

template <typename A1, typename A2, typename A3>
struct G4GPRLinearHierarchyT<G4GPRTypeList_3(A1, A2, A3)> : public A1, public A2, public A3 {};

template <typename A1, typename A2, typename A3, typename A4>
struct G4GPRLinearHierarchyT<G4GPRTypeList_4(A1, A2, A3, A4)> : public A1, public A2, public A3, public A4 {};

template <typename A1, typename A2, typename A3, typename A4, typename A5>
struct G4GPRLinearHierarchyT<G4GPRTypeList_5(A1, A2, A3, A4, A5)> : public A1, public A2, public A3, public A4, public A5 {};

template <typename A1, typename A2, typename A3,
	  typename A4, typename A5, typename A6>
struct G4GPRLinearHierarchyT<G4GPRTypeList_6(A1, A2, A3, A4, A5, A6)> : public A1, public A2, public A3, public A4, public A5, public A6 {};

template <typename A1, typename A2, typename A3,
	  typename A4, typename A5, typename A6, typename A7>
struct G4GPRLinearHierarchyT<G4GPRTypeList_7(A1, A2, A3, A4, A5, A6, A7)> : public A1, public A2, public A3, public A4, public A5, public A6, public A7 {};


template <typename A1, typename A2, typename A3,
	  typename A4, typename A5, typename A6, typename A7, typename A8>
struct G4GPRLinearHierarchyT<G4GPRTypeList_8(A1, A2, A3, A4, A5, A6, A7, A8)> : public A1, public A2, public A3, public A4, public A5, public A6, public A7, public A8 {};


#endif
