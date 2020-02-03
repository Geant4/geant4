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
// G4VDivisionParameterisation
//
// Class description:
//
// Base class for parameterisations defining divisions of volumes
// for different kind of CSG and specific solids.

// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// 21.04.10 - M.Asai, Added gaps
//---------------------------------------------------------------------
#ifndef G4VDIVISIONPARAMETERISATION_HH
#define G4VDIVISIONPARAMETERISATION_HH 1

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4VPVParameterisation.hh"
#include "G4RotationMatrix.hh"

enum DivisionType { DivNDIVandWIDTH, DivNDIV, DivWIDTH };

class G4VPhysicalVolume;
class G4VSolid;

class G4VDivisionParameterisation : public G4VPVParameterisation
{ 
  public:  // with description
  
    G4VDivisionParameterisation( EAxis axis, G4int nDiv, G4double width,
                                 G4double offset, DivisionType divType,
                                 G4VSolid* motherSolid = nullptr);
    virtual ~G4VDivisionParameterisation();
  
    virtual G4VSolid* ComputeSolid(const G4int, G4VPhysicalVolume*);

    virtual void ComputeTransformation(const G4int copyNo,
                                       G4VPhysicalVolume *physVol) const = 0;
  
    inline const G4String& GetType() const;
    inline EAxis GetAxis() const;
    inline G4int GetNoDiv() const;
    inline G4double GetWidth() const;
    inline G4double GetOffset() const;
    inline G4VSolid* GetMotherSolid() const;
    inline void SetType(const G4String& type);
    inline G4int VolumeFirstCopyNo() const;
    inline void SetHalfGap(G4double hg);
    inline G4double GetHalfGap() const;

  protected:  // with description

    void ChangeRotMatrix( G4VPhysicalVolume* physVol,
                          G4double rotZ = 0.0 ) const;
  
    G4int CalculateNDiv( G4double motherDim, G4double width,
                         G4double offset ) const;
    G4double CalculateWidth( G4double motherDim, G4int nDiv,
                             G4double offset ) const;
  
    virtual void CheckParametersValidity();
    void CheckOffset( G4double maxPar );
    void CheckNDivAndWidth( G4double maxPar );
    virtual G4double GetMaxParameter() const = 0;
    G4double OffsetZ() const;

  protected:
  
    G4String ftype;
    EAxis faxis;
    G4int fnDiv = 0;
    G4double fwidth = 0.0;
    G4double foffset = 0.0;
    DivisionType fDivisionType;
    G4VSolid* fmotherSolid = nullptr;
    G4bool fReflectedSolid = false;
    G4bool fDeleteSolid = false;

    static G4ThreadLocal G4RotationMatrix* fRot;

    static const G4int verbose;
    G4int theVoluFirstCopyNo = 1;

    G4double kCarTolerance;

    G4double fhgap = 0.0;
};

#include "G4VDivisionParameterisation.icc"

#endif
