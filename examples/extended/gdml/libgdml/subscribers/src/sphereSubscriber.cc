#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "sphere.hh"

#include "G4VSolid.hh"
#include "G4Sphere.hh"

#include <iostream>

class sphereSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  sphereSubscriber() {
    Subscribe( "sphere" );
  }
  virtual ~sphereSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    std::cout << "SPHERE SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const sphere* obj = dynamic_cast<const sphere*>( object );    
        
        std::cout << "GOT SPHERE " << obj->get_name() << std::endl;
      
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
        sval = obj->get_startphi();
        sval += "*"+aunit;
        double startphi = calc->Eval( sval );
        sval = obj->get_deltaphi();
        sval += "*"+aunit;
        double deltaphi = calc->Eval( sval );
        sval = obj->get_starttheta();
        sval += "*"+aunit;
        double starttheta = calc->Eval( sval );
        sval = obj->get_deltatheta();
        sval += "*"+aunit;
        double deltatheta = calc->Eval( sval );
        
        std::cout << "rmin:       "  << obj->get_rmin()       << lunit
                  << " rmin: "       << rmin << std::endl;
        std::cout << "rmax:       "  << obj->get_rmax()       << lunit
                  << " rmax: "       << rmax << std::endl;
        std::cout << "startphi:   "  << obj->get_startphi()   << aunit
                  << " startphi: "   << startphi << std::endl;
        std::cout << "deltaphi:   "  << obj->get_deltaphi()   << aunit
                  << " deltaphi: "   << deltaphi << std::endl;
        std::cout << "starttheta: "  << obj->get_starttheta() << aunit
                  << " starttheta: " << starttheta << std::endl;
        std::cout << "deltatheta: "  << obj->get_deltatheta() << aunit
                  << " deltatheta: " << deltatheta << std::endl;

        G4VSolid* newobj = new G4Sphere( name, rmin,       rmax,
                                               startphi,   deltaphi,
                                               starttheta, deltatheta);
        
      
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

DECLARE_SUBSCRIBER_FACTORY(sphereSubscriber)

