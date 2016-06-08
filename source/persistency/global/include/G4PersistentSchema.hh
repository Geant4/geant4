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
// $Id: G4PersistentSchema.hh,v 1.5.8.1 2001/06/28 19:11:33 gunter Exp $
// GEANT4 tag $Name:  $
//
// Persistent-capable schema definitions required for ooddlx processor
//
// History:
// 19.11.99 Y.Morita - Created

#ifndef G4_PERSISTENTSCHEMA_HH
#define G4_PERSISTENTSCHEMA_HH

#ifdef OO_DDL_TRANSLATION

// make forward declaration of persistent classes
// in the name space "Geant4"
#pragma ooschema Geant4

class G4PDCofThisEvent;
class G4PVDigit;
// G4VDigi is embedded in G4PVDigit
class G4VDigi;
class G4PVDigitsCollection;
class G4PHCofThisEvent;
class G4PVHit;
// G4VHit is embedded in G4PVHit
class G4VHit;
class G4PVHitsCollection;
class G4PEvent;
class G4PPrimaryVertex;
class G4PPrimaryParticle;
class G4PAffineTransform;
class G4PVSolid;
class G4PGeometryObjectMap;
class G4PLogicalVolume;
class G4PPVParameterised;
class G4PPVPlacement;
class G4PPVReplica;
class G4PVPhysicalVolume;
class G4PBooleanSolid;
class G4PDisplacedSolid;
class G4PIntersectionSolid;
class G4PSubtractionSolid;
class G4PUnionSolid;
class G4PBox;
class G4PCSGSolid;
class G4PCons;
class G4PHype;
class G4PPara;
class G4PSphere;
class G4PTorus;
class G4PTrap;
class G4PTrd;
class G4PTubs;
class G4PRun;

// back to the default name space
#pragma ooschema

#endif

#endif /* G4_PERSISTENTSCHEMA_HH */
