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
// $Id: G4Field.hh,v 1.6 2001-11-08 17:31:07 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4Field
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.

// History:
// - Created:  John Apostolakis, 10.03.1997

#ifndef G4FIELD_HH
#define G4FIELD_HH

class G4Field
{
  public:  // with description

      virtual void  GetFieldValue( const  double Point[4],
					  double *Bfield ) const = 0;
      G4Field(){;}
      virtual ~G4Field(){;}

     // A field signature function that can be used to insure
     // that the Equation of motion object and the G4Field object
     // have the same "field signature"?
};


#endif /* G4FIELD_HH */







