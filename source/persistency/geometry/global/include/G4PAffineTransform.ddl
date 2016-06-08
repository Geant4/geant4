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
// $Id: G4PAffineTransform.ddl,v 1.4.4.1 2001/06/28 19:11:27 gunter Exp $
// GEANT4 tag $Name:  $
//
// persistenet class of G4AffineTransform
//
// History:
// 10.11.99 Y.Morita, initial creation

#ifndef G4PAFFINETRANSFORM_DDL
#define G4PAFFINETRANSFORM_DDL

#include "HepODBMS/odbms/HepODBMS.h"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

class G4AffineTransform;

class G4PAffineTransform
 : public HepPersObj
{
  public:
	G4PAffineTransform( const G4AffineTransform aTrans );
	~G4PAffineTransform();

        G4AffineTransform MakeTransientObject();

  private:
	G4Pdouble rxx,rxy,rxz;
	G4Pdouble ryx,ryy,ryz;
	G4Pdouble rzx,rzy,rzz;
	
	G4Pdouble tx,ty,tz;
};

#endif

