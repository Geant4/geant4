#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "isotope.hh"

#include "G4Isotope.hh"

#include <iostream>

class isotopeSubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const
  {
    return this;
  }

public:
  isotopeSubscriber() {
    Subscribe( "isotope" );
  }

  virtual ~isotopeSubscriber() {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "ISOTOPE SUBSCRIBER:: " << std::endl;
    
    const isotope* obj = 0;
    
    if( object != 0 ) {
      try {
        obj = dynamic_cast<const isotope*>(object);
        
        if( obj != 0 ) {
          //std::cout << "GOT ISOTOPE " << obj->get_name() << std::endl;

          GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
          // The big moment, we're gonna create an isotope
          double z = calc->Eval( obj->get_Z() );
          double n = calc->Eval( obj->get_N() );
          std::string sA = obj->get_atom().get_value();
          sA += "*";
          sA += obj->get_atom().get_unit();
          double a = calc->Eval( sA );
          //std::cout << "Z: " << obj->get_Z() << " " << z << std::endl;
          //std::cout << "N: " << obj->get_N() << " " << n << std::endl;
          //std::cout << "A: " << sA << " " << a/g/mole << std::endl;
          G4Isotope* inew = new G4Isotope( obj->get_name(), z, n, a );
          std::cout << *inew << std::endl;
        }
      } catch(...) {
        std::cerr << "GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "GOT ZERO DATA POINTER!" << std::endl;
    }
    delete obj;
  }
  
private:
  double     m_z;
  double     m_n;
  double     m_a;
  // Geant4 does not require material properties to be specified for isotopes 
};

DECLARE_SUBSCRIBER_FACTORY(isotopeSubscriber)
