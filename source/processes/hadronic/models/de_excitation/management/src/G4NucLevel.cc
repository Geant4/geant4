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
// $Id: G4NucLevel.cc 88375 2015-02-16 17:31:21Z vnivanch $
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

G4NucLevel::G4NucLevel(size_t ntrans, G4double tgamma,
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
    for(size_t i=0; i<length; ++i) {
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
  for(size_t i=0; i<length; ++i) {
    delete fShellProbability[i];
  }
}

#ifdef G4VERBOSE
void G4NucLevel::PrintError(size_t idx, const G4String& ss) const
{
  G4cout << "G4NucLevel::PrintError: length= " << length << G4endl;
  for(size_t i=0; i<length; ++i) {
    G4cout << i << ". " /*<< fFinalIndex[i] << "  "*/ << fTrans[i]
	   << fGammaCumProbability[i] << " " << fTimeGamma << " "
           << fGammaProbability[i] << " "
	   << fMpRatio[i] << G4endl;
  }
  G4String sss = "G4NucLevel::"+ss+"()";
  G4ExceptionDescription ed;
  ed << "Index of a level " << idx << " >= " 
     << length << " (number of transitions)";
  G4Exception(sss,"had061",JustWarning,ed,"");
  throw G4HadronicException(__FILE__, __LINE__,"FATAL Hadronic Exception");
}  
#endif
