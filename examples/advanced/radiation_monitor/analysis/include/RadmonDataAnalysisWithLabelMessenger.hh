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
// File name:     RadmonDataAnalysisWithLabelMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelMessenger.hh,v 1.2 2006-06-28 13:44:09 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONDATAANALYSISWITHLABELMESSENGER_HH
 #define  RADMONDATAANALYSISWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonDataAnalysisWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonDataAnalysisWithLabelMessenger * Instance(void);
  
   void                                         AddAvailableDataAnalysis(const G4String & name);
   void                                         RemoveAvailableDataAnalysis(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonDataAnalysisWithLabelMessenger();
                                                RadmonDataAnalysisWithLabelMessenger(const RadmonDataAnalysisWithLabelMessenger & copy);
                                               ~RadmonDataAnalysisWithLabelMessenger();
   RadmonDataAnalysisWithLabelMessenger &       operator=(const RadmonDataAnalysisWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailableDataAnalyses;
   
  // Private variables
   AvailableDataAnalyses                        availableDataAnalyses;
   
   static RadmonDataAnalysisWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONDATAANALYSISWITHLABELMESSENGER_HH */
