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
// File name:     RadmonApplicationAnalysisSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationAnalysisSetup.hh,v 1.2 2006-06-28 13:44:47 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application analysis setup
//

#ifndef   RADMONAPPLICATIONANALYSISSETUP_HH
 #define  RADMONAPPLICATIONANALYSISSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonDataAnalysisWithLabelFactory;

 class RadmonApplicationAnalysisSetup
 {
  public:
   inline                                       RadmonApplicationAnalysisSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationAnalysisSetup();

  protected:
   G4bool                                       CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationAnalysisSetup();
                                                RadmonApplicationAnalysisSetup(const RadmonApplicationAnalysisSetup & copy);
   RadmonApplicationAnalysisSetup &             operator=(const RadmonApplicationAnalysisSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationAnalysisSetup.icc"
#endif /* RADMONAPPLICATIONANALYSISSETUP_HH */
