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
/// \file hadronic/Hadr02/include/G4GlaubAADataSetHandler.hh
/// \brief Definition of the G4GlaubAADataSetHandler class
//
#ifndef G4GlaubAADataSetHandler_h
#define G4GlaubAADataSetHandler_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4GlaubAADataSetHandler.hh
//
// Version:             0.A
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4Isotope.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4GlauberDataSet;
class G4DPMJET2_5Model;

#include "G4GlaubAADataSet.hh"
#include "G4FullGlaubAADataSet.hh"
#include "G4ParamType1GlaubAADataSet.hh"
#include "G4DPMJET2_5Model.hh"

typedef std::map< G4int, G4ParamType1GlaubAADataSet *, std::less< G4int > >
  G4GlaubAADataSetIndex;
////////////////////////////////////////////////////////////////////////////////
//
//class G4GlaubAADataSetHandler : public G4VGlauberDataSetHandler
class G4GlaubAADataSetHandler
{
  
  protected:
    G4GlaubAADataSetHandler ();
    
  public:
    virtual ~G4GlaubAADataSetHandler ();
    static G4GlaubAADataSetHandler* getInstance ();

    G4bool IsGlauberDataSetAvailable (const G4int AP, const G4int AT) const;
    G4bool SetCurrentGlauberDataSet (const G4int AP, const G4int AT,
      const G4double ppn = 0.0);
    G4GlaubAADataSet *GetCurrentGlauberDataSet () const;
    void   ResetCurrentGlauberDataSet ();
    G4double GetValueN (const G4double v, const G4double ppn1 = 0.0);
    G4double GetValueM (const G4double v, const G4double ppn1 = 0.0);

    void   SetMaxGlauberDataSets (const G4int n);
    G4int  GetMaxGlauberDataSets () const;
    
    void SetVerboseLevel (const G4int i);
    G4int GetVerboseLevel () const;
    G4int UnloadAllGlauberData ();
    
  private:
    G4VGlauberDataSet *LoadGlauberDataReturnPtr(const G4int AP, const G4int AT);

    G4bool   CheckIfSpace () const;
    G4int    GetIndexID (const G4int AP, const G4int AT) const;
    G4String GetStringID (const G4int AP, const G4int AT) const;
    G4String GetProjectileStringID (const G4int AP) const;
    G4String GetTargetStringID (const G4int AT) const;
    
    static G4GlaubAADataSetHandler *instance;

    G4GlaubAADataSetIndex          theIndex;
    G4String                       glauberDataSetDir;

    G4int                          maxGlauberDataSets;
    G4int                          cntGlauberDataSets;

    G4GlaubAADataSet              *theCurrentGlauberDataSet;
    G4double                       ppnCurrent;
    G4bool                         usingLocalGlauberDataSet;

    G4int                          maxArray;
    
    G4int                          verboseLevel;
};

inline void G4GlaubAADataSetHandler::SetVerboseLevel (const G4int i)
  {verboseLevel = i;}
  
inline G4int G4GlaubAADataSetHandler::GetVerboseLevel () const
  {return verboseLevel;}
#endif
