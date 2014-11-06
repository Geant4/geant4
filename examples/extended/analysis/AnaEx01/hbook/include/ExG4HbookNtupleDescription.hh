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
// $Id: G4HnInformation.hh 66310 2012-12-17 11:56:35Z ihrivnac $

/// \file hbook/include/ExG4HbookNtupleDescription.hh
/// \brief Definition of the ExG4HbookNtupleDescription structure

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookNtupleDescription_h
#define ExG4HbookNtupleDescription_h 1

#include "tools/hbook/wntuple"
#include "tools/ntuple_booking"

#include <map>

/// Structure containing the information related to one Hbook ntuple

struct ExG4HbookNtupleDescription
{
  ExG4HbookNtupleDescription() 
    :  fNtuple(0),
       fNtupleBooking(0),
       fNtupleIColumnMap(),
       fNtupleFColumnMap(),
       fNtupleDColumnMap() {}
 
  ~ExG4HbookNtupleDescription() 
     { delete fNtupleBooking; } 

  tools::hbook::wntuple* fNtuple; 
  tools::ntuple_booking*  fNtupleBooking; 
  std::map<G4int, tools::hbook::wntuple::column<int>* >    fNtupleIColumnMap;           
  std::map<G4int, tools::hbook::wntuple::column<float>* >  fNtupleFColumnMap;           
  std::map<G4int, tools::hbook::wntuple::column<double>* > fNtupleDColumnMap;           
};

#endif  

#endif  
