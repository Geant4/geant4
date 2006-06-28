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
//
//
// File name:     RadmonVDataAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysis.hh,v 1.3 2006-06-28 13:44:15 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
