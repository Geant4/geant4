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
/*
 * G4Serialize.hh
 *
 *  Created on: Jul 8, 2015
 *      Author: mkaramit
 */

#ifndef SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MOLECULES_MANAGEMENT_INCLUDE_G4SERIALIZE_HH_
#define SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MOLECULES_MANAGEMENT_INCLUDE_G4SERIALIZE_HH_

#include "globals.hh"

//_____________________________________________________________________________

template<typename T>
void WRITE(std::ostream& out, const T& toBeSaved)
{
  out.write((char*)(&toBeSaved), sizeof(toBeSaved));
}

//_____________________________________________________________________________

template<typename T>
void READ(std::istream& in, T& toBeSaved)
{
  in.read((char*)(&toBeSaved), sizeof(toBeSaved));
}

//_____________________________________________________________________________

template<>
void WRITE<G4String>(std::ostream& out, const G4String& name);

template<>
void READ<G4String>(std::istream& in, G4String& name);

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MOLECULES_MANAGEMENT_INCLUDE_G4SERIALIZE_HH_ */
