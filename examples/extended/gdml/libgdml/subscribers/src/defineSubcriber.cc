#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "define.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include <iostream>

class defineSubscriber : virtual public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const  {
    return this;
  }

public:
  defineSubscriber()  {
    Subscribe( "define" );
  }

  virtual ~defineSubscriber()  {
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object )  {
    //std::cout << "define SUBSCRIBER:: " << std::endl;
    GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
    const define* define_element = 0;
    if( object != 0 )    {
      try {
        define_element = dynamic_cast<const define*>(object);

        if( define_element->size() > 0 ) {
          // There is something defined in there
          unsigned int num = define_element->size();
          // Now let's traverse the define's structure
          for( unsigned int idx = 0; idx < num; idx++ ) {
            if( define_element->get_content(idx)->type() == SAXObject::contentGroup ) {
              // Fine that's what's expected here
              ContentChoice* choice = dynamic_cast<ContentChoice*>(define_element->get_content(idx));
              ContentGroup::ContentItem item = choice->content();
              if( item.tag == "constant" ) {
                define::constant* c = dynamic_cast<define::constant*>(item.object);
                //std::cout << "GOT constant " << c->get_name() << " = " << c->get_value() << std::endl;
                calc->RegisterConstant( c );
              } else if( item.tag == "quantity" ) {
                define::quantity* q = dynamic_cast<define::quantity*>(item.object);
                //std::cout << "GOT quantity " << q->get_name() << " = " << q->get_value() << q->get_unit() << std::endl;
                calc->RegisterPhysConstant( q );
              } else if( item.tag == "expression" ) {
                define::expression* e = dynamic_cast<define::expression*>(item.object);
                calc->RegisterExpression( e );
                //std::cout << "GOT expression " << e->get_name()
                //          << " = "             << e->get_text()
                //          << std::endl;
              } else if( item.tag == "position" ) {
                define::position* p = dynamic_cast<define::position*>(item.object);
                //std::cout << "GOT position " << p->get_name() << " = ("
                //                             << p->get_x()    << ","
                //                             << p->get_y()    << ","
                //                             << p->get_z()    << ")"
                //                             << p->get_unit() << std::endl;
                std::string
                expr = p->get_x() + "*" + p->get_unit();
                double dx = calc->Eval( expr.c_str() );
                expr = p->get_y() + "*" + p->get_unit();
                double dy = calc->Eval( expr.c_str() );
                expr = p->get_z() + "*" + p->get_unit();
                double dz = calc->Eval( expr.c_str() );
                G4ThreeVector* posobj = new G4ThreeVector( dx, dy, dz );
                GDMLProcessor::GetInstance()->AddPosition( p->get_name(), posobj );
              } else if( item.tag == "rotation" ) {
                define::rotation* r = dynamic_cast<define::rotation*>(item.object);
                //std::cout << "GOT rotation " << r->get_name() << " = " << "("
                //                             << r->get_x()    << ","
                //                             << r->get_y()    << ","
                //                             << r->get_z()    << ")"
                //                             << r->get_unit() << std::endl;
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
//                 HepTransform3D* result = new HepRotateX3D(dx);
//                 *result = *result * HepRotateY3D(dy);
//                 *result = *result * HepRotateZ3D(dz);
//                 G4RotationMatrix rm4 = result->getRotation();
//                 std::cout << "GOT rotation " << "rm4" << " = " << "("
//                                              << rm4.getPhi() << ","
//                                              << rm4.getTheta() << ","
//                                              << rm4.getPsi() << ")"
//                                              << r->get_unit() << std::endl;
              } else {
                // Problem, nothing else can appear in define element
                std::cerr << "define SUBSCRIBER:: Incorrect content model found!\a" << std::endl;
              }
            } else if( define_element->get_content(idx)->type() == SAXObject::element ) {
              // Whoops, this should not happen as the only content of define element is 0 or more choices
              std::cerr << "define SUBSCRIBER:: Incorrect content model found!\a" << std::endl;
            }
          }
        }
      } catch(...) {
        std::cerr << "GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
      delete object;
    }
  }
};

DECLARE_SUBSCRIBER_FACTORY(defineSubscriber)

