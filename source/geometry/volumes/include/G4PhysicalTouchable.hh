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
// $Id: G4PhysicalTouchable.hh,v 1.2 2005-03-03 17:09:12 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4PhysicalTouchable
//
// Class description:
//
// This is a concrete class that combines a simple physical volume,  
// with a touchable for the parent volume
//
// Behaviour expected: 
//    - pass all 'VPhysicalVolume' methods to pCurrentVol
//    - respond to GetParentTouchable method
//
// History:
// 14.02.05 J.Apostolakis Created

#ifndef G4PHYSICALTOUCHABLE_HH
#define G4PHYSICALTOUCHABLE_HH

#include "G4Types.hh"
#include "G4String.hh"

#include "geomdefs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4LogicalVolume;
class G4VPVParameterisation;
class G4VTouchable;
class G4NavigationHistory; 

class G4PhysicalTouchable : public G4VPhysicalVolume
{
  public:  // with description

    G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, 
			 const G4VTouchable* pParentTouchable); 
    G4PhysicalTouchable( G4VPhysicalVolume* pCurrentVol, 
			 const G4NavigationHistory& parentHist); 
      // Constructors, with existing touchable or history for use in one.

    virtual ~G4PhysicalTouchable();
      // Destructor, will be subclassed. Removes volume from volume Store.

    inline const G4VTouchable* GetParentTouchable() const; 
    void SetParentTouchable( const G4VTouchable* newParentT );  
      // Provide / set parent touchable
    inline G4VPhysicalVolume* GetCurrentVolume(); 
    inline const G4VPhysicalVolume* GetCurrentVolume() const; 
      // Access 

    inline void SetCurrentVolume( G4VPhysicalVolume* pCurrentVol ); 
      // Revise current volume pointer

    // inline G4VPhysicalVolume* operator ->(); 
      // Refer other methods to physical volume ???

    G4int GetMultiplicity() const { return fpPhysVol->GetMultiplicity(); } 
    G4bool IsMany() const { return fpPhysVol->IsMany(); } 

    G4int GetCopyNo() const { return fpPhysVol->GetCopyNo(); } 
    void  SetCopyNo(G4int CopyNo) { fpPhysVol->SetCopyNo(CopyNo); } 
      // Get/Set the volumes copy number.
    G4bool IsReplicated() const { return fpPhysVol->IsReplicated(); } 
      // Return true if replicated (single object instance represents
      // many real volumes), else false.
    G4bool IsParameterised() const { return fpPhysVol->IsParameterised(); }
      // Return true if parameterised (single object instance represents
      // many real parameterised volumes), else false.
    G4VPVParameterisation* GetParameterisation() const 
      { return fpPhysVol->GetParameterisation(); }
      // Return replicas parameterisation object 
    void GetReplicationData(EAxis& axis,
                                    G4int& nReplicas,
                                    G4double& width,
                                    G4double& offset,
                                    G4bool& consuming) const 
      { fpPhysVol->GetReplicationData(axis, nReplicas, width, offset, consuming); } 
      // Return replication information. No-op for no replicated volumes.

  private: 
    void CopyAttributes( G4VPhysicalVolume* pCurrentVol );
    G4PhysicalTouchable( const G4PhysicalTouchable & ) ; 
  
  private:
    G4VPhysicalVolume*   fpPhysVol; 
      // Current volume pointer

    const G4VTouchable*  fpTouchable;
    G4bool fCreatedParentTouch; 

}; 

inline G4VPhysicalVolume* G4PhysicalTouchable::GetCurrentVolume()
{
    return fpPhysVol;
}
inline const G4VPhysicalVolume* G4PhysicalTouchable::GetCurrentVolume() const
{
    return fpPhysVol;
}

inline const G4VTouchable* G4PhysicalTouchable::GetParentTouchable() const
      // Provide touchable
{
   return fpTouchable; 
}

#if 0
// Extra Inline methods
inline G4VPhysicalVolume* G4PhysicalTouchable::operator ->() // const
{
    return( fpPhysVol ? fpPhysVol : 0 );
}
#endif

  
#endif
