#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "setup.hh"

#include <iostream>

class setupSubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const
  {
    return this;
  }

public:
  setupSubscriber() {
    Subscribe( "setup" );
  }

  virtual ~setupSubscriber() {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object )
  {
    //std::cout << "SETUP SUBSCRIBER:: " << std::endl;
    
    const setup* obj = 0;
    
    if( object != 0 )
    {
      try {
        obj = dynamic_cast<const setup*>(object);       
        if( obj != 0 ) {
          //std::cout << "GOT SETUP " << obj->get_name()
          //          << " version "  << obj->get_version() << std::endl;
          
          if( obj->get_match() ) {
            // The selected setup has been found
            // Check & retrieve world volume
            setup::world w = obj->get_world();
            G4LogicalVolume* wlv = const_cast<G4LogicalVolume*>
                                           (GDMLProcessor::GetInstance()->
                                                           GetLogicalVolume(w.get_ref()));
            if( wlv == 0 ) {              
              std::cerr << "VOLUME SUBSCRIBER:: world volume "
                        << w.get_ref() << " not found!" << std::endl;
              std::cerr << "World volume " << obj->get_name() << " can't be created!" << std::endl;
              std::cerr << "Please, re-order your volumes or add the missing one..."
                        << std::endl;
              G4Exception( "Shutting-down due to error(s) in GDML input..." );                
            }
    		    G4VPhysicalVolume* wpv = new G4PVPlacement( 0, G4ThreeVector(),
                                                       wlv, wlv->GetName(),
                                                       0, false, 0 );
            GDMLProcessor::GetInstance()->SetWorldVolume( wpv );	
          }
        }
      } catch(...) {
        std::cerr << "GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "GOT ZERO DATA POINTER!" << std::endl;
    }
    delete obj;
  } 
};

DECLARE_SUBSCRIBER_FACTORY(setupSubscriber)
