#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "tube.hh"

#include "G4VSolid.hh"
#include "G4Tubs.hh"

#include <iostream>

class tubeSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  tubeSubscriber() {
    Subscribe( "tube" );
  }
  virtual ~tubeSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "TUBE SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const tube* obj = dynamic_cast<const tube*>( object );    
        
        //std::cout << "GOT TUBE " << obj->get_name() << std::endl;
      
        GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
      
        std::string lunit = obj->get_lunit();
        std::string aunit = obj->get_aunit();
        const std::string& name = obj->get_name();
        
        std::string sval = obj->get_rmin();
        sval += "*"+lunit;
        double rmin = calc->Eval( sval );
        sval = obj->get_rmax();
        sval += "*"+lunit;
        double rmax = calc->Eval( sval );
        sval = obj->get_z();
        sval += "*"+lunit;
        double dz = calc->Eval( sval ); dz = dz/2.;
        sval = obj->get_startphi();
        sval += "*"+aunit;
        double startphi = calc->Eval( sval );
        sval = obj->get_deltaphi();
        sval += "*"+aunit;
        double deltaphi = calc->Eval( sval );
        
//         std::cout << "rmin:       "  << obj->get_rmin()       << lunit
//                   << " rmin:      "  << rmin << std::endl;
//         std::cout << "rmax:       "  << obj->get_rmax()       << lunit
//                   << " rmax:      "  << rmax << std::endl;
//         std::cout << "z:          "  << obj->get_z()       << lunit
//                   << " dz:        "  << dz << std::endl;
//         std::cout << "startphi:   "  << obj->get_startphi()   << aunit
//                   << " startphi:  "  << startphi << std::endl;
//         std::cout << "deltaphi:   "  << obj->get_deltaphi()   << aunit
//                   << " deltaphi:  "  << deltaphi << std::endl;

        G4VSolid* newobj = new G4Tubs( name, rmin, rmax, dz, startphi, deltaphi );
      
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

DECLARE_SUBSCRIBER_FACTORY(tubeSubscriber)

