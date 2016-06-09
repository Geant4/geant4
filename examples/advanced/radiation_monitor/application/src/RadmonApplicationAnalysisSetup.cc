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
// File name:     RadmonApplicationAnalysisSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationAnalysisSetup.cc,v 1.3.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-02 $
//

// Include files
#include "RadmonApplicationAnalysisSetup.hh"
#include "RadmonApplicationOptions.hh"

#ifdef    G4ANALYSIS_USE
 #include "RadmonDataAnalysisDepositedEnergy.hh"

 #include "RadmonDataAnalysisWithLabelFactory.hh"


 #define DECLARE_ANALYSIS_CONSTRUCTOR(name)     constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendDataAnalysisWithLabel(constructor)

 G4bool RadmonApplicationAnalysisSetup :: CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * factory)
 {
  RadmonVDataAnalysisWithLabel * constructor;

  DECLARE_ANALYSIS_CONSTRUCTOR(RadmonDataAnalysisDepositedEnergy);

  return true;
 }
#else  /* G4ANALYSIS_USE */
 G4bool RadmonApplicationAnalysisSetup :: CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * /* factory */)
 {
  return false;
 }
#endif /* G4ANALYSIS_USE */
