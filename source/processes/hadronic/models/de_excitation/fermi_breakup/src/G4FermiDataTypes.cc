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
//
// G4FermiBreakUpAN alternative FermiBreakUp model
// by A. Novikov (January 2025)
//

#include "G4FermiDataTypes.hh"

std::string std::to_string(G4FermiAtomicMass mass)
{
  return std::to_string(G4FermiAtomicMass::ValueType(mass));
}

std::string std::to_string(G4FermiChargeNumber charge)
{
  return std::to_string(G4FermiChargeNumber::ValueType(charge));
}

std::ostream& std::operator<<(std::ostream& out, const G4FermiAtomicMass& mass)
{
  out << G4FermiAtomicMass::ValueType(mass);
  return out;
}

std::istream& std::operator>>(std::istream& in, G4FermiAtomicMass& mass)
{
  G4FermiAtomicMass::ValueType val;
  in >> val;
  mass = G4FermiAtomicMass(val);
  return in;
}

std::ostream& std::operator<<(std::ostream& out, const G4FermiChargeNumber& charge)
{
  out << G4FermiChargeNumber::ValueType(charge);
  return out;
}

std::istream& std::operator>>(std::istream& in, G4FermiChargeNumber& charge)
{
  G4FermiChargeNumber::ValueType val;
  in >> val;
  charge = G4FermiChargeNumber(val);
  return in;
}
