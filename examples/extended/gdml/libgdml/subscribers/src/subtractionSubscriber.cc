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
// $Id: subtractionSubscriber.cc,v 1.2 2002-06-03 12:09:35 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "BooleanSolidTypeSubscriber.hh"

#include "boolean_subtraction.hh"

#include "G4SubtractionSolid.hh"

class subtractionSubscriber : public BooleanSolidTypeSubscriber
{
public:
  subtractionSubscriber() {
    Subscribe( "subtraction" );
  }
  virtual ~subtractionSubscriber() {
  }
  
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    
    const boolean_subtraction* bs = dynamic_cast<const boolean_subtraction*>( object );
    
    // First, let the base class extract the data
    BooleanSolidTypeSubscriber::Activate( object );
    
    // At this moment both constituent solids are ready
    // Let's check whether we need to apply any transformations
    G4VSolid* solid_subtraction = 0;
    bool useTransform = false;
    G4ThreeVector translat;
    G4RotationMatrix transrot;
    
    if( m_translation != 0 ) {
      translat = *m_translation;
      useTransform = true;
    }
    
    if( m_rotation != 0 ) {
      transrot = *m_rotation;
      useTransform = true;
    }

    if( useTransform ) {    
      G4Transform3D transform( transrot, translat );
      solid_subtraction = new G4SubtractionSolid( bs->get_name(), m_first, m_second, transform );
    } else {      
      solid_subtraction = new G4SubtractionSolid( bs->get_name(), m_first, m_second );
    }
    
    GDMLProcessor::GetInstance()->AddSolid( bs->get_name(), solid_subtraction );
    delete object;
  }
};

DECLARE_SUBSCRIBER_FACTORY(subtractionSubscriber)
