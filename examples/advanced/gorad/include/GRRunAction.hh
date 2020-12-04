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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRRunAction.hh
//   Header file of Gorad Run Action class that takes care of
//   handling histograms and n-tuple.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRRunAction_h
#define GRRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class GRRunActionMessenger;
class GRRunAction;
#include "GRRun.hh"

#include <map>
class G4VPrimitivePlotter;

class GRHistoType
{
  friend class GRRunAction;
  friend class GRRun;
  private:
    GRHistoType()
    {;}
  private:
    G4int histID = -1;
    G4int histType = -1;
    G4int histDup = 1;

    G4int collID = -1;
    G4String meshName = "dummy";
    G4String primName = "dummy";
    G4int idx = -1;

    G4int collID2 = -1;
    G4String meshName2 = "dummy";
    G4String primName2 = "dummy";
    G4int idx2 = -1;

    G4int biasf = 0;
    G4double fuct = 1.;
    G4VPrimitivePlotter* pplotter = nullptr;
};
    
class GRRunAction : public G4UserRunAction
{
  friend class GRRun;
  public:
    GRRunAction();
    virtual ~GRRunAction();

    virtual G4Run* GenerateRun()
    { return new GRRun(this); }
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  private:
    GRRunActionMessenger* messenger;

  public:
    void SetVerbose(G4int);
    void ListHistograms();
    G4bool Open(G4int);
    G4bool SetAllPlotting(G4bool val=true);
    G4bool SetPlotting(G4int,G4bool val=true);
    void Flush();
    void Reset();

    G4int Create1D(G4String&,G4String&,G4int);
    G4int Create1DForPrimary(G4String&,G4bool);
    G4int Create1DForPlotter(G4String&,G4String&,G4bool);
    G4bool Set1D(G4int,G4int,G4double,G4double,G4String&,G4String&,G4bool);
    G4bool Set1DTitle(G4int,G4String&,G4String&,G4String&);
    G4bool Set1DYAxisLog(G4int,G4bool);

    G4int Create1P(G4String&,G4String&,G4int);
    G4bool Set1P(G4int,G4double,G4double,G4String&,G4String&,G4String&,G4String&);
    G4bool Set1PTitle(G4int,G4String&,G4String&,G4String&);

    G4int NtupleColumn(G4String&,G4String&,G4String&,G4int);

  private:
    void OpenFile();
    void DefineNTColumn();
    void MergeNtuple();

  public:
    inline void SetFileName(G4String& fn)
    { fileName = fn; }
    inline const G4String& GetFileName() const
    { return fileName; }
    inline G4int GetVerbose() const
    { return verbose; }
    inline void SetCarry(G4bool val = true)
    { ifCarry = val; }
    inline G4bool GetCarry() const
    { return ifCarry; }
    inline void SetOffset(G4int offset,G4int factor)
    { 
      id_offset = offset;
      id_factor = factor;
    }
    inline void GetOffset(G4int& offset,G4int& factor) const
    {
      offset = id_offset;
      factor = id_factor;
    }

  private:
    G4String fileName;
    G4bool fileOpen;
    G4int verbose;
    G4bool ifCarry;
    G4int id_offset;
    G4int id_factor;

  private:
    std::map<G4int,GRHistoType*> IDMap;
    std::map<G4int,GRHistoType*> NTMap;
};

#endif
