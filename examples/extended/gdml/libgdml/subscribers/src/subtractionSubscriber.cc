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
