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
// G4ReflectionFactory
//
// Class description:
//
// Class providing functions for volumes placements with a general
// transfomation that can contain a reflection.
// Reflection is then applied to a solid: a new G4ReflectedSolid
// instance is created and is placed with a transformation containing
// pure rotation and translation only.
// The pair of constituent and reflected logical volumes is
// considered as a generalised logical volume that is addressed
// by user specifying the constituent logical volume.

// Author: Ivana Hrivnacova (IN2P3/IJCLab Orsay), 16.10.2001
// --------------------------------------------------------------------
#ifndef G4_REFLECTION_FACTORY_HH
#define G4_REFLECTION_FACTORY_HH

#include "G4Types.hh"
#include "G4Transform3D.hh"
#include "geomdefs.hh"

#include <map>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4VPVDivisionFactory;

using G4PhysicalVolumesPair = std::pair<G4VPhysicalVolume*, G4VPhysicalVolume*>;
using G4ReflectedVolumesMap = std::map<G4LogicalVolume*, G4LogicalVolume*,  
                                       std::less<G4LogicalVolume*> >;
/**
 * @brief G4ReflectionFactory provides functions for volumes placements with
 * a general transfomation that can contain a reflection.
 * Reflection is then applied to a solid: a new G4ReflectedSolid instance is
 * created and is placed with a transformation containing pure rotation and
 * translation only.
 * The pair of constituent and reflected logical volumes is considered as a
 * generalised logical volume that is addressed by user specifying the
 * constituent logical volume.
 *
 * Decomposition of a general transformation that can include reflection
 * in a "reflection-free" transformation:
 * 
 * x(inM') = TG*x(inM)         TG - general transformation
 *         = T*(R*x(inM))      T  - "reflection-free" transformation
 *         = T* x(inReflM)   
 *
 * Daughters transformation:
 * When a volume V containing daughter D with transformation TD
 * is placed in mother M with a general tranformation TGV,
 * the TGV is decomposed. New reflected volume ReflV containing
 * a new daughter ReflD with reflected transformation ReflTD is created:
 * 
 * x(inV) = TD * x(inD);
 * x(inM) = TGV * x(inV) 
 *        = TV * R * x(inV) 
 *        = TV * R * TD * x(inD)
 *        = TV * R*TD*R-1 * R*x(inD)
 *        = TV * ReflTD * x(inReflD)
 */

class G4ReflectionFactory 
{
  using LogicalVolumesMapIterator = G4ReflectedVolumesMap::const_iterator;

  public:
  
    /**
     * Destructor.
     */
    ~G4ReflectionFactory();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4ReflectionFactory(const G4ReflectionFactory&) = delete;
    G4ReflectionFactory& operator=(const G4ReflectionFactory&) = delete;
 
    /**
     * Gets pointer to the instance of the singleton.
     */
    static G4ReflectionFactory* Instance();

    /**
     * Evaluates the passed transformation; if it contains reflection
     * it performs its decomposition, creates new reflected solid and
     * logical volume (or retrieves them from a map if the reflected
     * objects were already created), transforms the daughters (if present)
     * and places it in the given mother.
     *  @param[in] transform3D The transformation in 3D space that can contain
     *             a reflection.
     *  @param[in] pName The volume name.
     *  @param[in] LV Pointer to the logical volume to be placed.
     *  @param[in] motherLV Pointer to the logical volume of the mother.
     *  @param[in] isMany Not used.
     *  @param[in] copyNo The optional custom copy number.
     *  @param[in] surfCheck Boolean flag, if true activates check for overlaps
     *             with existing volumes (false by default).
     *  @returns A pair of physical volumes; the second physical volume
     *           is a placement in a reflected mother or nullptr if mother
     *           logical volume was not reflected.
     */
    G4PhysicalVolumesPair Place(const G4Transform3D& transform3D,
                                const G4String&  name,
                                      G4LogicalVolume* LV,
                                      G4LogicalVolume* motherLV,
                                      G4bool isMany, 
                                      G4int  copyNo,
                                      G4bool surfCheck = false);

    /**
     * Creates a replica in the given mother.
     *  @param[in] name The volume name.
     *  @param[in] LV Pointer to the logical volume to be replicated.
     *  @param[in] motherLV Pointer to the logical volume of the mother.
     *  @param[in] axis The axis along which do the replication.
     *  @param[in] nofReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the replicated object along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns A pair of physical volumes; the second physical volume is
     *           a replica in a reflected mother or nullptr if the mother
     *           logical volume was not reflected.
     */
    G4PhysicalVolumesPair Replicate(const G4String& name, 
                                          G4LogicalVolume* LV,
                                          G4LogicalVolume* motherLV,
                                          EAxis axis, 
                                          G4int nofReplicas, 
                                          G4double width,
                                          G4double offset = 0.);

    /**
     * Methods to create a division in the given mother, along with the
     * possible specifications for creating a division.
     *  @returns A pair of physical volumes; the second physical volume is
     *           a division in a reflected mother or nullptr if mother
     *           logical volume was not reflected.
     */
    G4PhysicalVolumesPair Divide(const G4String& name, 
                                       G4LogicalVolume* LV,
                                       G4LogicalVolume* motherLV,
                                       EAxis axis, 
                                       G4int nofDivisions, 
                                       G4double width,
                                       G4double offset);
    G4PhysicalVolumesPair Divide(const G4String& name, 
                                       G4LogicalVolume* LV,
                                       G4LogicalVolume* motherLV,
                                       EAxis axis, 
                                       G4int nofDivisions, 
                                       G4double offset);
    G4PhysicalVolumesPair Divide(const G4String& name, 
                                       G4LogicalVolume* LV,
                                       G4LogicalVolume* motherLV,
                                       EAxis axis, 
                                       G4double width,
                                       G4double offset);

    /**
     * Verbosity control.
     */
    void SetVerboseLevel(G4int verboseLevel);
    G4int GetVerboseLevel() const;

    /**
     * Sets/returns the name extension for the reflected solids and
     * logical volumes.
     */
    void SetVolumesNameExtension(const G4String& nameExtension);
    const G4String& GetVolumesNameExtension() const;

    /**
     * Sets/gets the precision factor for the scale consistency check.
     * The default value is set to 10*kCarTolerance.
     */
    void SetScalePrecision(G4double scaleValue);
    G4double GetScalePrecision() const;

    /**
     * Returns the consituent volume of the given reflected volume.
     * Returns nullptr if the given reflected volume was not found.
     */
    G4LogicalVolume* GetConstituentLV(G4LogicalVolume* reflLV) const;

    /**
     * Returns the reflected volume of the given consituent volume.
     * Returns nullptr if the given volume was not reflected.
     */
    G4LogicalVolume* GetReflectedLV(G4LogicalVolume* lv) const;

    /**
     * Returns true if the given volume has been already reflected
     * (i.e. is in the map of constituent volumes).
     */
    G4bool IsConstituent(G4LogicalVolume* lv) const;

    /**
     * Returns true if the given volume is a reflected volume
     * (i.e. is in the map reflected  volumes).
     */
    G4bool IsReflected(G4LogicalVolume* lv) const;

    /**
     * Returns a handle to the internal map of volumes which have
     * been reflected, after that placement or replication is performed.
     */
    const G4ReflectedVolumesMap& GetReflectedVolumesMap() const;

    /**
     * Clears the maps of constituent and reflected volumes.
     *  @note To be used exclusively when volumes are removed from the stores.
     */
    void Clean();

  private:  

    /**
     * Private singleton constructor.
     */
    G4ReflectionFactory();

    /**
     * Gets/creates the reflected solid and logical volume
     * and copies + transforms logical volumes of daughters.
     */
    G4LogicalVolume* ReflectLV(G4LogicalVolume* LV, G4bool surfCheck = false);

    /**
     * Creates the reflected solid and logical volume
     * and adds the logical volumes pair in the maps.
     */
    G4LogicalVolume* CreateReflectedLV(G4LogicalVolume* LV);

    /**
     * Reflects daughters recursively.
     */
    void ReflectDaughters(G4LogicalVolume* LV,
                          G4LogicalVolume* refLV, G4bool surfCheck = false);

    /**
     * Copies and transforms daughter of G4PVPlacement type of
     * a constituent volume into a reflected volume.
     */
    void ReflectPVPlacement(G4VPhysicalVolume* PV,
                            G4LogicalVolume* refLV, G4bool surfCheck = false);

    /**
     * Copies and transforms daughter of G4PVReplica type of
     * a constituent volume into a reflected volume.
     */ 
    void ReflectPVReplica(G4VPhysicalVolume* PV, G4LogicalVolume* refLV);

    /**
     * Copies and transforms daughter of G4PVDivision type of
     * a constituent volume into a reflected volume.
     */
    void ReflectPVDivision(G4VPhysicalVolume* PV, G4LogicalVolume* refLV);

    /**
     * Not yet implemented.
     * Should copy and transform daughter of G4PVParameterised type of
     * a constituent volume into a reflected volume. 
     */
    void ReflectPVParameterised(G4VPhysicalVolume* PV,
                                G4LogicalVolume* refLV, G4bool surfChk = false);

    /**
     * Returns true if the scale is negative, false otherwise.
     */
    G4bool IsReflection(const G4Scale3D& scale) const;

    /**
     * Checks if scale correspond to fScale, if not gives exception.
     */
    void CheckScale(const G4Scale3D& scale) const;

    /**
     * Checks if the division factory is instanciated, if not gives exception.
     */
    G4VPVDivisionFactory* GetPVDivisionFactory() const;

    /**
     * Dump method, for debugging purpose.
     */
    void PrintConstituentLVMap() const;  
  
  private:

    static G4ThreadLocal G4ReflectionFactory* fInstance;
    static const G4String fDefaultNameExtension;
    static const G4Scale3D fScale;
    G4double fScalePrecision;

    G4int fVerboseLevel = 0;
    G4String fNameExtension;
    G4ReflectedVolumesMap fConstituentLVMap;
    G4ReflectedVolumesMap fReflectedLVMap;
};

#endif
