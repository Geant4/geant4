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
// File: G4VPDigitsCollectionIO.cc
//
// History:
//   '01.08.16  Youhei Morita  Initial creation

#include "G4VPDigitsCollectionIO.hh"

// Implementation of Constructor #1
G4VPDigitsCollectionIO::G4VPDigitsCollectionIO( G4std::string detName,
                                                G4std::string colName )
 : m_verbose(0), f_detName(detName), f_colName(colName)
{}

// Implementation of operator== 
G4bool G4VPDigitsCollectionIO::operator== (const G4VPDigitsCollectionIO& right) const
{
  return ( (f_detName == right.f_detName) &&
           (f_colName == right.f_colName) );
}

// End of G4VPDigitsCollectionIO.cc

