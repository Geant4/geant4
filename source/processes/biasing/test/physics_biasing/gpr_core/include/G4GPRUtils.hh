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
// $Id: G4GPRUtils.hh,v 1.4 2007-09-06 22:10:10 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRUTILS_HH
#define G4GPRUTILS_HH

#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRTypeList.hh"

namespace G4GPRUtils {

  template <typename Result, typename A1>
  void Operator(Result* result, G4GPRLinearHierarchyT<G4GPRTypeList_1(A1)>* input)
  {
    input->A1::operator()(result);
  }
  
  template <typename Result, typename A1, typename A2>
  void Operator(Result* result, G4GPRLinearHierarchyT<G4GPRTypeList_2(A1, A2)>* input)
  {
    input->A1::operator()(result);
    input->A2::operator()(result);
  }

  template <typename Result, typename A1, typename A2, typename A3>
  void Operator(Result* result, G4GPRLinearHierarchyT<G4GPRTypeList_3(A1, A2, A3)>* input)
  {
    input->A1::operator()(result);
    input->A2::operator()(result);
    input->A3::operator()(result);
  }

  template <typename Result, typename A1, typename A2, typename A3, typename A4>
  void Operator(Result* result, G4GPRLinearHierarchyT<G4GPRTypeList_4(A1, A2, A3, A4)>* input)
  {
    input->A1::operator()(result);
    input->A2::operator()(result);
    input->A3::operator()(result);
    input->A4::operator()(result);
  }
}

#endif
