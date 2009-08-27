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
// Description:
//
// A file for link definitions of ROOT. More information here:
// http://root.cern.ch/drupal/content/selecting-dictionary-entries-linkdefh
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class G4TModelParams+;
#pragma link C++ class G4TSimulationTool+;
#pragma link C++ class G4TSimHelper+;
#pragma link C++ class G4TAnalysisTool+;
#pragma link C++ class G4TPlotHelper+;
#pragma link C++ class G4TData+;
#pragma link C++ class G4TDataItem+;
#pragma link C++ class G4TParticlesDAL+;
#pragma link C++ class G4TCatalog+;
#pragma link C++ class G4TTool+;
#pragma link C++ class G4TDataBase+;
#pragma link C++ class G4TGUIHelper+;
#pragma link C++ class G4TSimulationGUI+;
#pragma link C++ class G4TAnalysisGUI+;
#pragma link C++ class G4TPublicationGUI+;
#pragma link C++ class G4TToolMenuGUI+;



#pragma link C++ global gPlotHelper;
#pragma link C++ global gParticlesDAL;
#pragma link C++ global gCatalog;
#pragma link C++ global gSimHelper;
#pragma link C++ global gSimulationTool;
#pragma link C++ global gAnalysisTool;
#pragma link C++ global gTestDB;
#pragma link C++ global gGUIHelper;

#endif



