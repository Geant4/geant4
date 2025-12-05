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
// G4VPVDivisionFactory
//
// Class description:
//
// Abstract factory that defines interfaces for creating volume divisions. 

// Author: Ivana Hrivnacova (IJCLab, Orsay), 04.05.2004
// ------------------------------------------------------------------------
#ifndef G4VPVDIVISIONFACTORY_HH
#define G4VPVDIVISIONFACTORY_HH

#include "geomdefs.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VPVParameterisation;

/**
 * @brief G4VPVDivisionFactory is an abstract factory that defines the
 * interfaces for creating volume divisions.
 */

class G4VPVDivisionFactory
{
  public:

    /**
     * Default Destructor.
     */
    virtual ~G4VPVDivisionFactory() = default;
    
    /**
     * Returns a pointer to the unique instance of the factory.
     */
    static G4VPVDivisionFactory* Instance();
    
    /**
     * Creates a division, using number of divisions and width.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns The pointer to physical divided volume.
     */
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double width,
                             const G4double offset ) = 0;

    /**
     * Creates a division, using the number of divisions.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns The pointer to physical divided volume.
     */
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const EAxis pAxis,
                             const G4int nReplicas,
                             const G4double offset ) = 0;

    /**
     * Creates a division, using the width of the division slice.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns The pointer to physical divided volume.
     */
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const EAxis pAxis,
                             const G4double width,
                             const G4double offset ) = 0;

    /**
     * Creates a division, using a parameterisation algorithm.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] param The pointer to the parameterisation algorithm.
     *  @returns The pointer to physical divided volume.
     */
    virtual G4VPhysicalVolume* CreatePVDivision(
                             const G4String& pName,
                                   G4LogicalVolume* pLogical,
                                   G4LogicalVolume* pMother,
                             const G4VPVParameterisation* param) = 0;

    /**
     * Returns true if 'pv' is a pointer to a division.
     */
    virtual G4bool IsPVDivision(const G4VPhysicalVolume* pv) const = 0;

  protected:

    /**
     * Protected default Constructor.
     */
    G4VPVDivisionFactory() = default;

  protected:

    /** The static pointer of the singleton instance. */
    static G4ThreadLocal G4VPVDivisionFactory* fgInstance;
};

#endif
