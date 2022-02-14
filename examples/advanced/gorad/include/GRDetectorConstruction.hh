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
// GRDetectorConstruction.hh
//   Header file of the detector construction. It reads a GDML file.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRDetectorConstruction_H
#define GRDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4GDMLParser;
class GRDetectorConstructionMessenger;
class GRGeomImpBiasWorld;

class GRDetectorConstruction : public G4VUserDetectorConstruction
{
  friend class GRGeomImpBiasWorld;

  public:
    GRDetectorConstruction();
    virtual ~GRDetectorConstruction();
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDAndField();

  private:
    GRDetectorConstructionMessenger* messenger;
    G4GDMLParser* parser;
    G4String gdmlFile;
    G4VPhysicalVolume* fWorld;
    G4bool initialized;

    static G4double worldSize;

  public:
    G4bool SetGDMLFile(G4String&);
    const G4String& GetGDMLFile() 
    { return gdmlFile; }

    static G4double GetWorldSize()
    { return worldSize; }

  private:
    void Read();

  public:
    void ListSolids(G4int);
    void ListLogVols(G4int);
    void ListPhysVols(G4int);
    void ListRegions(G4int);
    G4bool CreateRegion(G4String&,G4String&);
    G4bool CheckOverlap(G4String&,G4int,G4int,G4double);

  public:
    void ListAllMaterial();
    G4bool ListMaterial(G4String&);
    void DumpNistMaterials();
    G4bool CreateMaterial(G4String&);
    G4bool GetMaterial(G4String&);
    G4int SetMaterial(G4String&,G4String&);

  private:
    G4bool applyGeomImpBias = false;
    GRGeomImpBiasWorld* geomImpBiasWorld = nullptr;
    struct GeomImpParameters
    {
      G4double radius = -1.;
      G4ThreeVector pos0;
      G4int nLayer = 0;
      G4double radiusT = -1.;
      G4ThreeVector posT;
      G4int factor = 2;
      G4double prob = 1.;
    } geoImpP;
      
  public:
    void GeomImp(G4int n,G4double r)
    {
      applyGeomImpBias = true;
      geoImpP.nLayer = n;
      geoImpP.radius = r;
    }
    G4bool ApplyGeomImpBias() const
    { return applyGeomImpBias; }
    void GeomImpLocate(G4ThreeVector p0)
    { geoImpP.pos0 = p0; }
    G4double GeomImpInnerRadius(G4double rt)
    {
      if(rt>geoImpP.radius)
      { return -rt; }
      geoImpP.radiusT = rt;
      return rt;
    }
    G4double GeomImpLocateTgt(G4ThreeVector pT)
    {
      G4ThreeVector dp = geoImpP.pos0 - pT;
      G4double rt = geoImpP.radiusT;
      if(rt<0.)
      { rt = geoImpP.radius / geoImpP.nLayer; }
      if(dp.mag() >= (geoImpP.radius - rt))
      { return (geoImpP.radius - rt); }
      geoImpP.posT = pT;
      return 0.;
    }
    void GeomImpFactor(G4int f)
    { geoImpP.factor = f; }
    void GeomImpProb(G4double p)
    { geoImpP.prob = p; }
};

#endif
