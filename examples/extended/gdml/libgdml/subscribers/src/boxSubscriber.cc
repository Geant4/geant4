#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "box.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"

#include <iostream>

class boxSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  boxSubscriber() {
    Subscribe( "box" );
  }
  virtual ~boxSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "BOX SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const box* obj = dynamic_cast<const box*>( object );    
        
        //std::cout << "GOT BOX " << obj->get_name() << std::endl;
      
        GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
      
        const std::string& name = obj->get_name();
        std::string sval = obj->get_x();
        sval += "*";
        sval += obj->get_lunit();
        double dx = calc->Eval( sval ); dx = dx/2.;
        sval = obj->get_y();
        sval += "*";
        sval += obj->get_lunit();
        double dy = calc->Eval( sval ); dy = dy/2.;
        sval = obj->get_z();
        sval += "*";
        sval += obj->get_lunit();
        double dz = calc->Eval( sval ); dz = dz/2.;
        
        //std::cout << "x: " << obj->get_x() << obj->get_lunit() << " dx: " << dx << std::endl;
        //std::cout << "y: " << obj->get_y() << obj->get_lunit() << " dy: " << dy << std::endl;
        //std::cout << "z: " << obj->get_z() << obj->get_lunit() << " dx: " << dz << std::endl;

        G4VSolid* newobj = new G4Box( name, dx, dy, dz );
      
        GDMLProcessor::GetInstance()->AddSolid( name, newobj );      
      } catch(...) {
        std::cerr << "GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "GOT ZERO DATA POINTER!" << std::endl;
    }
    delete object;
  }
};

DECLARE_SUBSCRIBER_FACTORY(boxSubscriber)

