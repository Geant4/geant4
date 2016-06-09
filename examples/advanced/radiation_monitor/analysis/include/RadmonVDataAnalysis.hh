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
// File name:     RadmonVDataAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysis.hh,v 1.4 2006/06/29 16:07:15 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Abstract class of an analysis piece
//

#ifndef   RADMONVDATAANALYSIS_HH
 #define  RADMONVDATAANALYSIS_HH

 // Include files
 #include "RadmonSensitiveDetectorDataStorer.hh"
 #include "RadmonHit.hh"
 #include "G4String.hh"
 #include "AIDA/ITuple.h" 
 
 class RadmonVDataAnalysis : public RadmonSensitiveDetectorDataStorer
 {
  public:
   inline virtual                              ~RadmonVDataAnalysis();
    
   virtual void                                 SetDataAnalysisAttribute(const G4String & attributeName, const G4String &value) = 0;

   virtual G4String                             ObtainColumnsDeclaration(const G4String & prefix) = 0;
   virtual void                                 InitializeFromTuple(const G4String & prefix, const AIDA::ITuple * tuple) = 0;
   virtual void                                 StoreIntoTuple(RadmonHitsCollection * hitsCollection, AIDA::ITuple * tuple) = 0;
   
  protected:
   inline                                       RadmonVDataAnalysis();

  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysis(const RadmonVDataAnalysis & copy);
   RadmonVDataAnalysis &                        operator=(const RadmonVDataAnalysis & copy);
 };
 
 // Inline implementations
 #include "RadmonVDataAnalysis.icc"
#endif /* RADMONVDATAANALYSIS_HH */
