#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "trd.hh"

#include "G4VSolid.hh"
#include "G4Trd.hh"

#include <iostream>

class trdSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  trdSubscriber() {
    Subscribe( "trd" );
  }
  virtual ~trdSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "TRAPEZOID SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const trd* obj = dynamic_cast<const trd*>( object );    
        
        //std::cout << "GOT TRAPEZOID " << obj->get_name() << std::endl;
      
        GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
      
        const std::string& name = obj->get_name();
        std::string lunit = obj->get_lunit();
        
        std::string sval = obj->get_x1();
        sval += "*"+lunit;
        double dx1 = calc->Eval( sval ); dx1 = dx1/2.;
        sval = obj->get_y1();
        sval += "*"+lunit;
        double dy1 = calc->Eval( sval ); dy1 = dy1/2.;
        sval = obj->get_x2();
        sval += "*"+lunit;
        double dx2 = calc->Eval( sval ); dx2 = dx2/2.;
        sval = obj->get_y2();
        sval += "*"+lunit;
        double dy2 = calc->Eval( sval ); dy2 = dy2/2.;
        sval = obj->get_z();
        sval += "*"+lunit;
        double dz = calc->Eval( sval ); dz = dz/2.;
        
//         std::cout << "x1: " << obj->get_x1() << lunit << " dx1: " << dx1 << std::endl;
//         std::cout << "y1: " << obj->get_y1() << lunit << " dy1: " << dy1 << std::endl;
//         std::cout << "x2: " << obj->get_x2() << lunit << " dx2: " << dx2 << std::endl;
//         std::cout << "y2: " << obj->get_y2() << lunit << " dy2: " << dy2 << std::endl;
//         std::cout << "z:  " << obj->get_z()  << lunit << " dx:  " << dz  << std::endl;

        G4VSolid* newobj = new G4Trd( name, dx1, dx2, dy1, dy2, dz);
        
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

DECLARE_SUBSCRIBER_FACTORY(trdSubscriber)

