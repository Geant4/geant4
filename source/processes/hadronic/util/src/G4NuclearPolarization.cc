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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      File name:     G4NuclearPolarization
//
//      Author:        Jason Detwiler (jasondet@gmail.com)
// 
//      Creation date: Aug 2015
//
//      V.Ivanchenko left only polarization tensor and add operators
//
// -------------------------------------------------------------------

#include "G4NuclearPolarization.hh"
#include <iomanip>

G4NuclearPolarization::G4NuclearPolarization()
{
  Unpolarize();
}

G4NuclearPolarization::~G4NuclearPolarization()
{}

G4bool G4NuclearPolarization::operator==(const G4NuclearPolarization &right) const
{
  return (fPolarization == right.fPolarization);
}

G4bool G4NuclearPolarization::operator!=(const G4NuclearPolarization &right) const
{
  return (fPolarization != right.fPolarization);
}

std::ostream& operator<<(std::ostream& out, const G4NuclearPolarization& p)
{
  out << " P = [ {";
  size_t kk = p.fPolarization.size();
  for(size_t k=0; k<kk; ++k) {
    if(k>0) { out << "       {"; }
    size_t kpmax = (p.fPolarization[k]).size();
    for(size_t kappa=0; kappa<kpmax; ++kappa) {
      if(kappa > 0) { out << "}  {"; }
      out << p.fPolarization[k][kappa].real() << " + " 
	  << p.fPolarization[k][kappa].imag() << "*i";
    }
    if(k+1 < kk) { out << "}" << G4endl; }
  }
  out << "} ]" << G4endl;
  return out; 
}
