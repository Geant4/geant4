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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4GlaubAADataSetHandler.cc
/// \brief Implementation of the G4GlaubAADataSetHandler class
//
// $Id: G4GlaubAADataSetHandler.cc 77519 2013-11-25 10:54:57Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4GlaubAADataSetHandler.cc
//
// Version:             0.A
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET


#include "G4GlaubAADataSetHandler.hh"
#include "G4FullGlaubAADataSet.hh"
#include "G4ParamType1GlaubAADataSet.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4StableIsotopes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4HadronicException.hh"

#include "G4DPMJET2_5Interface.hh"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

G4GlaubAADataSetHandler* G4GlaubAADataSetHandler::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
G4GlaubAADataSetHandler::G4GlaubAADataSetHandler ()
{
//
//
// theIndex is a std::map of the Glauber data loaded.  It uses a unique
// identifier based on the projectile and target A & Z to point to each Glauber
// data set object.
//
  theIndex.clear();
//
//
// By default, there is no limit on the number of datasets to load.
//
  maxGlauberDataSets =-1;
  cntGlauberDataSets = 0;
//
//
// The glauber data sets are assumed to be in Geant4 data directory for
// DPMJET II.5, defined by an environment variable.
//
  if ( !getenv("G4DPMJET2_5DATA") )
  {
    G4cout <<"ENVIRONMENT VARIABLE G4DPMJET2_5DATA NOT SET " <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
      "Please setenv G4DPMJET2_5DATA to point to the dpmjet2.5 data files.");
  }
  glauberDataSetDir = getenv("G4DPMJET2_5DATA");
//
//
// Initialise the local copy of the Glauber data.
//
  theCurrentGlauberDataSet = 0;
  ppnCurrent               = 0.0;
  usingLocalGlauberDataSet = false;
  
  maxArray   = 200;
//
//
// No verbose output by default.
//
  SetVerboseLevel(0);
}
////////////////////////////////////////////////////////////////////////////////
//
G4GlaubAADataSetHandler::~G4GlaubAADataSetHandler ()
{
  UnloadAllGlauberData ();
}
////////////////////////////////////////////////////////////////////////////////
//
G4GlaubAADataSetHandler* G4GlaubAADataSetHandler::getInstance ()
{
//
//
// THERE CAN BE ONLY ONE!! .... I mean, one G4GlaubAADataSetHandler.
//
  if (instance == 0) instance = new G4GlaubAADataSetHandler();
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int 
G4GlaubAADataSetHandler::GetIndexID (const G4int AP, const G4int AT) const
{
  return AP*1000 + AT;
}
////////////////////////////////////////////////////////////////////////////////
//
G4String 
G4GlaubAADataSetHandler::GetProjectileStringID (const G4int AP) const
{
  std::ostringstream os;
  os <<"ap" <<AP;
  return os.str();
}
////////////////////////////////////////////////////////////////////////////////
//
G4String G4GlaubAADataSetHandler::GetTargetStringID (const G4int AT) const
{
  std::ostringstream os;
  os <<"at" <<AT;
  return os.str();
}

////////////////////////////////////////////////////////////////////////////////
//
G4String G4GlaubAADataSetHandler::GetStringID (const G4int AP, const G4int AT)
 const
{
  G4String ID = GetProjectileStringID(AP) + "_" + GetTargetStringID(AT);
  return ID;
}
////////////////////////////////////////////////////////////////////////////////
//
// LoadGlauberDataRtnPtr
//
// Function LoadGlauberData to load data for a specific projectile AP and target
// AT.  Note that, unlike in the previous implementation, there is no overwrite
// option since there is only one source of glaubAA data.
//
// Returns a pointer the the glauber data set object, if exists otherwise
// returns NULL.
//
G4VGlauberDataSet *G4GlaubAADataSetHandler::LoadGlauberDataReturnPtr
  (const G4int AP, const G4int AT)
{
  G4ParamType1GlaubAADataSet *glauberData = 0;
//
//
// Check whether an file with an appropriate name exists in the directory, 
// which will probably contain the Glauber data.
//
  G4String glauberID = GetStringID(AP,AT);
  G4String filename  = glauberDataSetDir + "/" + glauberID + ".glaubaadat";
  std::ifstream glauberFile(filename);
  if (glauberFile) {
//
//
// File exists, so read data into an appropriate object type.  Do a double-
// check and make sure we want to retain Glauber data relating to this target.
// Note that this GDS handler only deals with parameterised sets (type 1). Also
// note that peek returns an ASCII charater code, and therefore since values of
// '0' to '9' and 'A' to 'Z' are acceptable, these must be decoded.
//
    G4int i    = -1;
    G4int asci = glauberFile.peek();
    if (asci >= 48 && asci <= 57)      i = asci - 48;
    else if (asci >= 65 && asci <= 90) i = asci - 55;
    if (i!=1) {
      G4cerr <<"ERROR IN G4GlaubAADataSetHandler::LoadGlauberDataReturnPtr"
             <<G4endl;
      G4cerr <<"GLAUBER FILE " <<filename <<G4endl;
      G4cerr <<"IDENTIFIED AS GLAUBER FILE TYPE " <<i <<G4endl;
      G4cerr <<"THIS IS NOT SUPPORTED" <<G4endl;
      G4cerr <<G4endl;
    }
    else {
      glauberData = new G4ParamType1GlaubAADataSet();
      glauberFile >>(*glauberData);
//
//
// Here we check the A of the projectile and target...
//
      G4int AP1 = glauberData->GetAP();
      G4int AT1 = glauberData->GetAT();
      G4bool found = AP1==AP && AT1==AT;
      if (found) {
//
//
// Keep the object and record in the index, or if there is insufficient
// space, record it only locally and set the usingLocalGlauberDataSet
// flag.  In the latter case, the object is deleted after the call to
// G4DPMJET2_5Model.
//
        theCurrentGlauberDataSet = glauberData;
        if (CheckIfSpace()) {
          G4int n = GetIndexID(AP,AT);
          theIndex.insert (G4GlaubAADataSetIndex::value_type(n,glauberData));
          cntGlauberDataSets++;
          usingLocalGlauberDataSet = false;
          if (verboseLevel >=2) {
            G4cout <<"****************************************"
                   <<"****************************************"
                   <<G4endl;
            G4cout <<"In G4GlaubAADataSetHandler::LoadGlauberDataReturnPtr"
                   <<G4endl;
            G4cout <<"LOADED GLAUBER DATA SET FOR PROJECTILE A = " <<AP
                   <<" & TARGET A = " <<AT <<G4endl;
            G4cout <<"AS RECORD NUMBER = " <<theIndex.size() <<G4endl;
            G4cout <<"****************************************"
                   <<"****************************************"
                   <<G4endl;
          }
        }
        else {
          usingLocalGlauberDataSet = true;
          if (verboseLevel >=2) {
            G4cout <<"****************************************"
                   <<"****************************************"
                   <<G4endl;
            G4cout <<"In G4GlaubAADataSetHandler::LoadGlauberDataReturnPtr"
                   <<G4endl;
            G4cout <<"LOADED GLAUBER DATA SET FOR PROJECTILE A = " <<AP
                   <<" & TARGET A = " <<AT <<" TEMPORARILY" <<G4endl;
            G4cout <<"****************************************"
                   <<"****************************************"
                   <<G4endl;
          }
        }
      }
      else {
        G4cerr <<"WARNING IN G4GlaubAADataSetHandler::LoadGlauberData"
               <<G4endl;
        G4cerr <<"GLAUBER FILE " <<filename <<" LOOKED LIKE IT SHOULD CONTAIN"
               <<G4endl;
        G4cerr <<"DATA FOR AP = " <<AP <<" AND AT = " <<AT <<G4endl;
        G4cerr <<"BUT CONTAINED AP = " <<AP1 <<" AND AT = " <<AT1 <<G4endl;
        G4cerr <<G4endl;
        delete glauberData;
        glauberData = 0;
      }
    }
    glauberFile.close();
  }
    
  return glauberData;
}
////////////////////////////////////////////////////////////////////////////////
//
// UnloadAllGlauberData
//
// Member function to delete all objects containg Glauber data.  This is very
// easy since all of the objects are pointed to in the map "theIndex".
//
// The value returned is the number of data sets removed.
//
G4int G4GlaubAADataSetHandler::UnloadAllGlauberData ()
{
  G4int n = 0;

  for (G4GlaubAADataSetIndex::iterator it=theIndex.begin(); it!=theIndex.end();
    it++)
  {
    G4ParamType1GlaubAADataSet *glauberData = it->second;
    cntGlauberDataSets--;
    delete glauberData;
    n++;
  }
  theIndex.clear();

  return n;
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool G4GlaubAADataSetHandler::IsGlauberDataSetAvailable (const G4int AP,
  const G4int AT) const
{
  G4int glauberID                          = GetIndexID (AP,AT);
  G4GlaubAADataSetIndex::const_iterator it = theIndex.find(glauberID);
  if (it == theIndex.end()) {
    G4String glauberStrID = GetStringID(AP,AT);
    G4String filename     = glauberDataSetDir + "/" + glauberStrID +
                            ".glaubaadat";
    std::ifstream glauberFile;
    glauberFile.open(filename);
    if (glauberFile) {
      glauberFile.close();
      return true;
    }
    else {
      return false;
    }
  }
  
  return true;
}
////////////////////////////////////////////////////////////////////////////////
//
// CheckIfSpace
//
// This member function reads the contents of the file GluaberIndex.dat, 
// and then proceeds to read all Glauber functions. This can deal with both 
// full Glauber data sets and parameterised data sets.
//
G4bool G4GlaubAADataSetHandler::CheckIfSpace () const
{
  if ( maxGlauberDataSets >= 0 &&
       ( (G4int) theIndex.size() ) >= maxGlauberDataSets) {
    G4cout <<"WARNING: G4GlaubAADataSetHandler::CheckIfSpace:"
           <<G4endl;
    G4cout <<"MAXIMUM NUMBER OF GLAUBER DATASETS REACHED : "
           <<maxGlauberDataSets <<G4endl;;
    return false;
  } 

  return true;
}
////////////////////////////////////////////////////////////////////////////////
//
// SetMaxGlauberDataSets
//
// This member function controls how many data sets are loaded, in order to
// avoid overloading the memory with data.
//
void G4GlaubAADataSetHandler::SetMaxGlauberDataSets (const G4int n)
{
  if (n > -1 && cntGlauberDataSets > n) {
    G4cerr <<"ERROR IN G4GlaubAADataSetHandler::SetMaxGlauberDataSets"
           <<G4endl;
    G4cerr <<"MAXIMUM NUMBER OF GLAUBER DATA SETS WOULD BE EXCEEDED IF YOU"
           <<G4endl;
    G4cerr <<"SET THIS VALUE.  ATTEMPTED TO SET TO " <<n
           <<"WHEN VALUE ALREADY" <<cntGlauberDataSets
           <<G4endl;
    G4cerr <<G4endl;
  }
  else {
    maxGlauberDataSets = n;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
// SetCurrentGlauberDataSet
//
// Called from G4DPMJET2_5Model to identifies a pointer to the appropriate
// Glauber data set, based in the AP and AT provided. If the apprproate dataset
// is not already loaded, an attempt will be made to find it and load if there
// is space.
//
// true is returned if a data set has been found and loaded;
// false if returned otherwise.
//
G4bool G4GlaubAADataSetHandler::SetCurrentGlauberDataSet (const G4int AP,
  const G4int AT, const G4double ppn)
{
  G4int glauberID                    = GetIndexID (AP,AT);
  G4GlaubAADataSetIndex::iterator it = theIndex.find(glauberID);
  if (it == theIndex.end()) {
//
//
// Have we exceeded any limits for the maximum number of data sets loaded?
// If not, try to locate the data and load it.
//
    theCurrentGlauberDataSet = 
      (G4ParamType1GlaubAADataSet*) LoadGlauberDataReturnPtr(AP,AT);
    if (theCurrentGlauberDataSet != 0) {
      ppnCurrent        = ppn;

      dtumat_.rprojj[0] = theCurrentGlauberDataSet->rproj;
      dtumat_.rtagg[0]  = theCurrentGlauberDataSet->rtarg;
      dtumat_.bstepp[0] = theCurrentGlauberDataSet->bstep;
      dtumat_.bmaxx[0]  = theCurrentGlauberDataSet->bmax;
      return true;
    }
    else {
      return false;
    }
  }
  else {
//
//
// The data are already loaded, so point to this.
//
    theCurrentGlauberDataSet = it->second;
    ppnCurrent        = ppn;

    dtumat_.rprojj[0] = theCurrentGlauberDataSet->rproj;
    dtumat_.rtagg[0]  = theCurrentGlauberDataSet->rtarg;
    dtumat_.bstepp[0] = theCurrentGlauberDataSet->bstep;
    dtumat_.bmaxx[0]  = theCurrentGlauberDataSet->bmax;
    return true;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4GlaubAADataSet *G4GlaubAADataSetHandler::GetCurrentGlauberDataSet () const
{
  return theCurrentGlauberDataSet;
}
////////////////////////////////////////////////////////////////////////////////
//
 void G4GlaubAADataSetHandler::ResetCurrentGlauberDataSet ()
{
  if (usingLocalGlauberDataSet) delete theCurrentGlauberDataSet;
  theCurrentGlauberDataSet = 0;
  ppnCurrent               = 0.0;
  usingLocalGlauberDataSet = false;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4GlaubAADataSetHandler::GetValueN (const G4double v,
  const G4double ppn1)
{
  G4double ppn = 0.0;
  if (ppn1*GeV > 1.0*eV)            ppn = ppn1;
  else if (ppnCurrent*GeV > 1.0*eV) ppn = ppnCurrent;
  else return 0.0;

  return 
    ((G4ParamType1GlaubAADataSet*) theCurrentGlauberDataSet)->GetValueN(v,ppn);
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4GlaubAADataSetHandler::GetValueM (const G4double v,
  const G4double ppn1)
{
  G4double ppn = 0.0;
  if (ppn1*GeV > 1.0*eV)            ppn = ppn1;
  else if (ppnCurrent*GeV > 1.0*eV) ppn = ppnCurrent;
  else return 0.0;

  return 
    ((G4ParamType1GlaubAADataSet*) theCurrentGlauberDataSet)->GetValueM(v,ppn);
}
////////////////////////////////////////////////////////////////////////////////
//
#endif
