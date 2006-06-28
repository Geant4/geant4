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
// File name:     RadmonDetectorMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisMessenger.hh,v 1.2 2006-06-28 13:43:51 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONANALYSISMESSENGER_HH
 #define  RADMONANALYSISMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVAnalysisLayout;
 
 class RadmonAnalysisMessenger : public RadmonMessenger
 {
  public:
                                                RadmonAnalysisMessenger(RadmonVAnalysisLayout * layout);
                                               ~RadmonAnalysisMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonAnalysisMessenger();
                                                RadmonAnalysisMessenger(const RadmonAnalysisMessenger & copy);
   RadmonAnalysisMessenger &                    operator=(const RadmonAnalysisMessenger & copy);

  // Private Attributes
   RadmonVAnalysisLayout *                      analysisLayout;

  // Commands
   RADMON_DECLARE_COMMAND(SetOutputFileName);
   RADMON_DECLARE_COMMAND(SetOutputFileFormat);
  
   RADMON_DECLARE_COMMAND(CreateSensitiveDetector);
   RADMON_DECLARE_COMMAND(SetSensitiveDetectorType);
   RADMON_DECLARE_COMMAND(RemoveSensitiveDetector);
   
   RADMON_DECLARE_COMMAND(CreateSensitiveDetectorType);
   RADMON_DECLARE_COMMAND(RemoveSensitiveDetectorType);

   RADMON_DECLARE_COMMAND(AppendDataAnalysisToSensitiveDetectorType);
   RADMON_DECLARE_COMMAND(SetDataAnalysisType);
   RADMON_DECLARE_COMMAND(RemoveDataAnalysis);

   RADMON_DECLARE_COMMAND(SetDataAnalysisAttribute);
   RADMON_DECLARE_COMMAND(ClearDataAnalysisAttribute);

   RADMON_DECLARE_COMMAND(DumpLayout);
   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONANALYSISMESSENGER_HH */
