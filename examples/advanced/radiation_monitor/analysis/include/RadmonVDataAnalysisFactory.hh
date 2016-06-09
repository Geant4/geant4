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
// File name:     RadmonVDataAnalysisFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysisFactory.hh,v 1.1.2.2 2006/06/29 16:07:19 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Abstract class of a factory of analysis items
//

#ifndef   RADMONVDATAANALYSISFACTORY_HH
 #define  RADMONVDATAANALYSISFACTORY_HH

 // Forward declaration
 class RadmonVDataAnalysis;
 class G4String;

 class RadmonVDataAnalysisFactory
 {
  public:
   inline virtual                              ~RadmonVDataAnalysisFactory();

   virtual RadmonVDataAnalysis *                CreateDataAnalysis(const G4String & dataAnalysis) = 0;

  protected:
   inline                                       RadmonVDataAnalysisFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysisFactory(const RadmonVDataAnalysisFactory & copy);
   RadmonVDataAnalysisFactory &                 operator=(const RadmonVDataAnalysisFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVDataAnalysisFactory.icc"
#endif /* RADMONVDATAANALYSISFACTORY_HH */
