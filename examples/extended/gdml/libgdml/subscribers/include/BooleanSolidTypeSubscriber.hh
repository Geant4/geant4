#ifndef BOOLEANSOLIDTYPESUBSCRIBER_H
#define BOOLEANSOLIDTYPESUBSCRIBER_H 1

#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "BooleanSolidType.hh"

#include "define.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4VSolid.hh"

#include <iostream>

class BooleanSolidTypeSubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  BooleanSolidTypeSubscriber()
    : m_translation( 0 ), m_rotation( 0 ), m_first( 0 ), m_second( 0 ) {
  }

  virtual ~BooleanSolidTypeSubscriber() {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "BOOLEAN SOLID SUBSCRIBER:: ";
    
    GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
    
    if( object != 0 ) {
      try {
        // Try the element type first
        const BooleanSolidType* obj = dynamic_cast<const BooleanSolidType*>( object );
        
        if( obj != 0 ) {
          //std::cout << "GOT BOOLEAN SOLID " << obj->get_name() << std::endl;

          const ContentSequence* bseq = obj->get_content();
          size_t seqlen               = bseq->size();
          
          for( size_t idx = 0; idx < seqlen; idx++ ) {
            // Loop over sequence of inner elements
            if( bseq->content(idx).tag == "first" ) {
              BooleanSolidType::first* s1 = dynamic_cast<BooleanSolidType::first*>
                                                        ( bseq->content(idx).object );
              std::string s1name = s1->get_ref();
              m_first = const_cast<G4VSolid*>( GDMLProcessor::GetInstance()->GetSolid( s1name ) );
              if( m_first == 0 ) {
                // We're in trouble we can't create boolean solid
                std::cerr << "BOOLEAN SOLID SUBSCRIBER:: solid "
                          << s1->get_ref() << " not found!" << std::endl;
                std::cerr << "Boolean solid " << obj->get_name() << " can't be created!" << std::endl;
                std::cerr << "Please, re-order your solids or add the missing one..."
                          << std::endl;
                G4Exception( "Shutting-down due to error(s) in GDML input..." );
              }
            } else if( bseq->content(idx).tag == "second" ) {
              BooleanSolidType::second* s2 = dynamic_cast<BooleanSolidType::second*>
                                                         ( bseq->content(idx).object );
              std::string s2name = s2->get_ref();
              m_second = const_cast<G4VSolid*>( GDMLProcessor::GetInstance()->GetSolid( s2name ) );
              if( m_second == 0 ) {
                // We're in trouble we can't create boolean solid
                std::cerr << "BOOLEAN SOLID SUBSCRIBER:: solid "
                          << s2->get_ref() << " not found!" << std::endl;
                std::cerr << "Boolean solid " << obj->get_name() << " can't be created!" << std::endl;
                std::cerr << "Please, re-order your solids or add the missing one..."
                          << std::endl;
                G4Exception( "Shutting-down due to error(s) in GDML input..." );
              }
            } else if( bseq->content(idx).tag == "choice" ) {
              // We apparently got either position or rotation or both
              const ContentChoice* bschoice = dynamic_cast<const ContentChoice*>
                                                          ( bseq->content(idx).object );
              if( bschoice->content().tag == "position" ) {
//                BooleanSolidType::position* p = dynamic_cast<BooleanSolidType::position*>
                define::position* p = dynamic_cast<define::position*>
                                                            ( bschoice->content().object );
//                 std::cout << "GOT position " << p->get_name() << " = ("
//                                              << p->get_x()    << ","
//                                              << p->get_y()    << ","
//                                              << p->get_z()    << ")"
//                                              << p->get_unit() << std::endl;
                std::string
                expr = p->get_x() + "*" + p->get_unit();
                double dx = calc->Eval( expr.c_str() );
                expr = p->get_y() + "*" + p->get_unit();
                double dy = calc->Eval( expr.c_str() );
                expr = p->get_z() + "*" + p->get_unit();
                double dz = calc->Eval( expr.c_str() );
                G4ThreeVector* posobj = new G4ThreeVector( dx, dy, dz );
                GDMLProcessor::GetInstance()->AddPosition( p->get_name(), posobj );
                m_translation = posobj;
              } else if( bschoice->content().tag == "rotation" ) {
//                BooleanSolidType::rotation* r = dynamic_cast<BooleanSolidType::rotation*>
                define::rotation* r = dynamic_cast<define::rotation*>
                                                            ( bschoice->content().object );
//                 std::cout << "GOT rotation " << r->get_name() << " = " << "("
//                                              << r->get_x() << ","
//                                              << r->get_y() << ","
//                                              << r->get_z() << ")"
//                                              << r->get_unit() << std::endl;
                std::string
                expr = r->get_x() + "*" + r->get_unit();
                double dx = calc->Eval( expr.c_str() );
                expr = r->get_y() + "*" + r->get_unit();
                double dy = calc->Eval( expr.c_str() );
                expr = r->get_z() + "*" + r->get_unit();
                double dz = calc->Eval( expr.c_str() );
                G4RotationMatrix* rotobj = new G4RotationMatrix;
                rotobj->rotateX( dx );
                rotobj->rotateY( dy );
                rotobj->rotateZ( dz );
                GDMLProcessor::GetInstance()->AddRotation( r->get_name(), rotobj );
                m_rotation = rotobj;
              } else if( bschoice->content().tag == "positionref" ) {
                BooleanSolidType::positionref* p = dynamic_cast<BooleanSolidType::positionref*>
                                                               ( bschoice->content().object );
                if( !p->get_ref().empty() ) {
			            m_translation = const_cast<G4ThreeVector*>( GDMLProcessor::GetInstance()->
                                                               GetPosition( p->get_ref() ) );
                } 
              } else if( bschoice->content().tag == "rotationref" ) {
                BooleanSolidType::rotationref* r = dynamic_cast<BooleanSolidType::rotationref*>
                                                               ( bschoice->content().object );
                if( !r->get_ref().empty() ) {
                  m_rotation = const_cast<G4RotationMatrix*>( GDMLProcessor::GetInstance()->
                                                               GetRotation( r->get_ref() ) );
                } 
              }
            }
          }
        }
      } catch(...) {
        std::cerr << "MATERIAL SUBSCRIBER:: GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "MATERIAL SUBSCRIBER:: GOT ZERO DATA POINTER!" << std::endl;
    }
  }
  
protected:
  G4ThreeVector*    m_translation;
  G4RotationMatrix* m_rotation;
  G4VSolid*         m_first;
  G4VSolid*         m_second;
};

#endif // BOOLEANSOLIDTYPESUBSCRIBER_H
