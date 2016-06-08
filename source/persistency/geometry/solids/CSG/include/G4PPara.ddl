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
// $Id: G4PPara.ddl,v 1.7 2001/07/11 10:02:23 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $

// Class Description:
//   Persistent version of G4PPara solid.

// History:
// 19.06.98 A.Kimura Converted G4Para.hh

#ifndef G4PPara_DDL
#define G4PPara_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4Para;
class G4VSolid;

#include "G4ThreeVector.hh"

class G4PPara
 : public G4PCSGSolid
{
public: // With description
    G4PPara(const G4Para* thePara);
    virtual ~G4PPara() ;
    // Constructor and Destructor

    G4VSolid* MakeTransientObject() const ;
    // Creates a transient boolean solid object.

public:
    G4ThreeVector GetSymAxis() const {
      G4double cosTheta
	    = 1.0/sqrt(1+fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi) ;
      return G4ThreeVector(fTthetaCphi*cosTheta,fTthetaSphi*cosTheta,cosTheta) ;
   }

public: // With description
    virtual G4GeometryType  GetEntityType() const {return G4String("G4Para");}
    // Returns the G4GeometryType of this solid.

private:
    G4double fDx,fDy,fDz;
    G4double fTalpha,fTthetaCphi,fTthetaSphi;
};
   	
#endif

