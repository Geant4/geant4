#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "volume.hh"
#include "child.hh"
#include "BooleanSolidType.hh"

#include "G4Material.hh"

#include "G4AssemblyVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include <iostream>

#include "g4std/strstream"

class volumeSubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  volumeSubscriber() {
    Subscribe( "volume" );
  }

  virtual ~volumeSubscriber() {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    std::cout << "VOLUME SUBSCRIBER:: " << std::endl;
    
    //GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
    GDMLProcessor* processor = GDMLProcessor::GetInstance();
    
    const volume* obj = 0;
            
    if( object != 0 ) {
      try {
        obj = dynamic_cast<const volume*>( object );
        
        if( obj != 0 ) {
          std::cout << "GOT VOLUME " << obj->get_name() << std::endl;

          // Let's analyze content model if volume
          const ContentSequence* seq = obj->get_content();
          size_t count = seq->size();                

          G4Material*      vmaterial = 0;
          G4VSolid*        vsolid    = 0;
          G4PVPlacement*   vchild    = 0;
          G4LogicalVolume* vnew      = 0;
          
		      for( size_t i = 0; i < count; i++ ) {
            if( seq->content(i).tag == "materialref" ) {
              // Check & retrieve material
              VolumeType::materialref* mref = dynamic_cast<VolumeType::materialref*>
                                                          ( seq->content(i).object );
              if( (vmaterial = G4Material::GetMaterial( mref->get_ref())) == 0 ) {
                std::cerr << "VOLUME SUBSCRIBER:: material "
                          << mref->get_ref() << " not found!" << std::endl;
                std::cerr << "Volume " << obj->get_name() << " can't be created!" << std::endl;
                std::cerr << "Please, re-order your materials or add the missing one..."
                          << std::endl;
                G4Exception( "Shutting-down due to error(s) in GDML input..." );                
              }
            } else if( seq->content(i).tag == "solidref" ) {
              // Check & retrieve solid
              VolumeType::solidref* sref = dynamic_cast<VolumeType::solidref*>
                                                          ( seq->content(i).object );
              vsolid = const_cast<G4VSolid*>( processor->GetSolid( sref->get_ref()) );
              if( vsolid == 0 ) {
                std::cerr << "VOLUME SUBSCRIBER:: solid "
                          << sref->get_ref() << " not found!" << std::endl;
                std::cerr << "Volume " << obj->get_name() << " can't be created!" << std::endl;
                std::cerr << "Please, re-order your solids or add the missing one..."
                          << std::endl;
                G4Exception( "Shutting-down due to error(s) in GDML input..." );                
              }
              // At this point we should have material & solid so we can create log. volume
              vnew = new G4LogicalVolume( vsolid, vmaterial, obj->get_name() );
       	      GDMLProcessor::GetInstance()->AddLogicalVolume( obj->get_name(), vnew );
            } else if( seq->content(i).tag == "child" ) {
              // Analyze each child's content
              child* c = dynamic_cast<child*>( seq->content(i).object );
              const ContentSequence* child_seq = c->get_content();
              size_t ccount = child_seq->size();
              SinglePlacementType::volumeref*   vr   = 0;
              BooleanSolidType::positionref*    pr   = 0;
              BooleanSolidType::rotationref*    rr   = 0;
              G4LogicalVolume*                  plog = 0;
              G4AssemblyVolume*                 alog = 0;
              G4ThreeVector*                    ppos = 0;
              G4RotationMatrix*                 prot = 0;
              bool doAssemblyInprint                 = false;
              for( size_t cidx = 0; cidx < ccount; cidx++ ) {
                if( child_seq->content(cidx).tag == "volumeref" ) {
                  // Check & retrieve volume
                  vr = dynamic_cast<SinglePlacementType::volumeref*>
                                   ( child_seq->content(cidx).object );
                  plog = const_cast<G4LogicalVolume*>
                                   (processor->GetLogicalVolume( vr->get_ref()));
                  if( plog == 0 ) {
                    // Let's check if an assembly request was ment
                    alog = const_cast<G4AssemblyVolume*>
                                     (processor->GetAssemblyVolume( vr->get_ref()));
                    if( alog == 0 ) {
                      std::cerr << "VOLUME SUBSCRIBER:: child volume "
                                << vr->get_ref() << " not found!" << std::endl;
                      std::cerr << "Volume " << obj->get_name() << " can't be created!" << std::endl;
                      std::cerr << "Please, re-order your volumes or add the missing one..."
                                << std::endl;
                      G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                    }
                    doAssemblyInprint = true;
                  }
                } else if( child_seq->content(cidx).tag == "positionref" ) {
                  // Check & retrieve position
                  pr = dynamic_cast<BooleanSolidType::positionref*>
                                   ( child_seq->content(cidx).object );
                  ppos = const_cast<G4ThreeVector*>
                                   ( processor->GetPosition( pr->get_ref() ) );
                  if( ppos == 0 ) {
                    std::cerr << "VOLUME SUBSCRIBER:: child's position "
                              << pr->get_ref() << " not found!" << std::endl;
                    std::cerr << "Volume " << obj->get_name() << " can't be created!" << std::endl;
                    std::cerr << "Please, check your position definitions or add the missing one..."
                              << std::endl;
                    G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                  }
                } else if( child_seq->content(cidx).tag == "rotationref" ) {
                  // Check & retrieve rotation
                  rr = dynamic_cast<BooleanSolidType::rotationref*>
                                   ( child_seq->content(cidx).object );
                  prot = const_cast<G4RotationMatrix*>
                                   ( processor->GetRotation( rr->get_ref() ) );
                  if( prot == 0 ) {
                    std::cerr << "VOLUME SUBSCRIBER:: child's rotation "
                              << rr->get_ref() << " not found!" << std::endl;
                    std::cerr << "Volume " << obj->get_name() << " can't be created!" << std::endl;
                    std::cerr << "Please, check your rotation definitions or add the missing one..."
                              << std::endl;
                    G4Exception( "Shutting-down due to error(s) in GDML input..." );                
                  }
                  // At this point we should have everything ready to create a child
                  if( doAssemblyInprint ) {
                    alog->MakeImprint( plog, *ppos, prot );
                  } else {
                    std::strstream pvname;
                    pvname << "pv_" << vr->get_ref() << "_" << (cidx-2) << std::ends;
                    vchild = new G4PVPlacement( prot, *ppos, plog, pvname.str(),
                                                             vnew, false, cidx-2 );
                    processor->AddPhysicalVolume( pvname.str(), vchild );	
                  }
                } else {
                  // Should not happen
                  ;
                }
              }
            } else {
              // Should not happen
              ;
            }
		      }              
        }
      } catch(...) {
        std::cerr << "VOLUME SUBSCRIBER:: GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "VOLUME SUBSCRIBER:: GOT ZERO DATA POINTER!" << std::endl;
    }
    delete object;
  }
};

DECLARE_SUBSCRIBER_FACTORY(volumeSubscriber)
