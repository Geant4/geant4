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
// $Id: assemblySubscriber.cc,v 1.3 2002-06-03 12:09:35 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "assembly.hh"
#include "child.hh"
#include "BooleanSolidType.hh"

#include "G4Material.hh"

#include "G4AssemblyVolume.hh"

#include <iostream>

#include "g4std/strstream"

class assemblySubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  assemblySubscriber() {
    Subscribe( "assembly" );
  }

  virtual ~assemblySubscriber() {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "ASSEMBLY VOLUME SUBSCRIBER:: " << std::endl;
    
    //GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
    GDMLProcessor* processor = GDMLProcessor::GetInstance();
    
    const assembly* obj = 0;
            
    if( object != 0 ) {
      try {
        obj = dynamic_cast<const assembly*>( object );
        
        if( obj != 0 ) {
          //std::cout << "GOT ASSEMBLY VOLUME " << obj->get_name() << std::endl;

          // Let's analyze content model if volume
          const ContentSequence* seq = obj->get_content();
          size_t count = seq->size();                

          G4AssemblyVolume* anew = new G4AssemblyVolume;
          
		      for( size_t i = 0; i < count; i++ ) {
            if( seq->content(i).tag == "child" ) {
              // Analyze each child's content
              child* c = dynamic_cast<child*>( seq->content(i).object );
              const ContentSequence* child_seq = c->get_content();
              size_t ccount = child_seq->size();
              SinglePlacementType::volumeref*   vr   = 0;
              BooleanSolidType::positionref*    pr   = 0;
              BooleanSolidType::rotationref*    rr   = 0;
              G4AssemblyTriplet*                ptrp = 0;
              G4LogicalVolume*                  plog = 0;
              G4ThreeVector*                    ppos = 0;
              G4RotationMatrix*                 prot = 0;
              for( size_t cidx = 0; cidx < ccount; cidx++ ) {
                if( child_seq->content(cidx).tag == "volumeref" ) {
                  // Check & retrieve volume
                  vr = dynamic_cast<SinglePlacementType::volumeref*>
                                   ( child_seq->content(cidx).object );
                  plog = const_cast<G4LogicalVolume*>(processor->GetLogicalVolume( vr->get_ref()));
                  if( plog == 0 ) {
                    std::cerr << "ASSEMBLY VOLUME SUBSCRIBER:: child volume "
                              << vr->get_ref() << " not found!" << std::endl;
                    std::cerr << "Assembly volume " << obj->get_name() << " can't be created!" << std::endl;
                    std::cerr << "Please, re-order your volumes or add the missing one..." << std::endl;
                    std::cerr << "\nNOTE! Assembly can't contain another assembly!\n" << std::endl;
                    G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                  }
                } else if( child_seq->content(cidx).tag == "positionref" ) {
                  // Check & retrieve position
                  pr = dynamic_cast<BooleanSolidType::positionref*>
                                   ( child_seq->content(cidx).object );
                  ppos = const_cast<G4ThreeVector*>( processor->GetPosition( pr->get_ref() ) );
                  if( ppos == 0 ) {
                    std::cerr << "ASSEMBLY VOLUME SUBSCRIBER:: child's position "
                              << pr->get_ref() << " not found!" << std::endl;
                    std::cerr << "Assembly volume " << obj->get_name() << " can't be created!" << std::endl;
                    std::cerr << "Please, check your position definitions or add the missing one..."
                              << std::endl;
                    G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                  }
                } else if( child_seq->content(cidx).tag == "rotationref" ) {
                  // Check & retrieve rotation
                  rr = dynamic_cast<BooleanSolidType::rotationref*>
                                   ( child_seq->content(cidx).object );
                  prot = const_cast<G4RotationMatrix*>( processor->GetRotation( rr->get_ref() ) );
                  if( prot == 0 ) {
                    std::cerr << "ASSEMBLY VOLUME SUBSCRIBER:: child's rotation "
                              << rr->get_ref() << " not found!" << std::endl;
                    std::cerr << "Assembly volume " << obj->get_name() << " can't be created!" << std::endl;
                    std::cerr << "Please, check your rotation definitions or add the missing one..."
                              << std::endl;
                    G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                  }
                  // At this point we should have everything ready to create assembly triplet
                  ptrp = new G4AssemblyTriplet( plog, *ppos, prot );
                  anew->AddPlacedVolume( plog, *ppos, prot );
                } else {
                  // Should not happen
                  ;
                }
              } // end of for(;;)
              processor->AddAssemblyVolume( obj->get_name(), anew );	
            } else {
              // Should not happen
              ;
            }
		      }
        }
      } catch(...) {
        std::cerr << "ASSEMBLY VOLUME SUBSCRIBER:: GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "ASSEMBLY VOLUME SUBSCRIBER:: GOT ZERO DATA POINTER!" << std::endl;
    }
    delete object;
  }
};

DECLARE_SUBSCRIBER_FACTORY(assemblySubscriber)
