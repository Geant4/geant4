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
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//  13.02.2015 Design change for gamma de-excitation 
//      
// -------------------------------------------------------------------

#include "G4NucLevel.hh"
#include "G4HadronicException.hh"
#include <iomanip>

G4NucLevel::G4NucLevel(std::size_t ntrans, G4double tgamma,
		       const std::vector<G4int>&   vTrans,
		       const std::vector<G4float>& wLevelGamma,
		       const std::vector<G4float>& wGamma,
      	               const std::vector<G4float>& vRatio,
		       const std::vector<const std::vector<G4float>*>& wShell)
  : length(ntrans), fTimeGamma(tgamma)
{
  if(0 < length) { 
    fTrans.reserve(length);
    fGammaCumProbability.reserve(length);
    fGammaProbability.reserve(length);
    fMpRatio.reserve(length);
    fShellProbability.reserve(length);
    for(std::size_t i=0; i<length; ++i) {
      fTrans.push_back(vTrans[i]);
      fGammaCumProbability.push_back(wLevelGamma[i]);
      fGammaProbability.push_back(wGamma[i]);
      fMpRatio.push_back(vRatio[i]);
      fShellProbability.push_back(wShell[i]);
    }
  }
}

G4NucLevel::~G4NucLevel()
{
  for(std::size_t i=0; i<length; ++i) {
    delete fShellProbability[i];
  }
}

void G4NucLevel::StreamInfo(std::ostream& out) const
{
  G4long prec = out.precision(4);
  for(std::size_t i=0; i<length; ++i) {
    out << std::setw(12) << FinalExcitationIndex(i) 
	<< std::setw(4) << TransitionType(i)
	<< std::setw(7) << fMpRatio[i] 
	<< std::setw(7) << fGammaCumProbability[i]
	<< std::setw(7) << fGammaProbability[i]
	<< "\n";
    const std::vector<G4float>* vec = fShellProbability[i];
    if(vec) {
      std::size_t len = vec->size();
      out << "              ";
      for(std::size_t j=0; j<len; ++j) { out << std::setw(7) << (*vec)[j]; }
      out << "\n";
    }
  }
  out.precision(prec);
}
