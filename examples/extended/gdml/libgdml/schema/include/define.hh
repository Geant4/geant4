#ifndef DEFINE_H
#define DEFINE_H 1


#include "SAXObject.hh"
#include "defineType.hh"
#include "IdentifiableConstantType.hh"
#include "IdentifiableQuantityType.hh"
#include "IdentifiableExpressionType.hh"
#include "positionType.hh"
#include "rotationType.hh"

class define : virtual public SAXObject, public defineType {
public:
  define() {
  }
  virtual ~define() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
public:
  class constant : public SAXObject, public IdentifiableConstantType {
    public:
      constant() {
      }
      virtual ~constant() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };
  class quantity : public SAXObject, public IdentifiableQuantityType {
    public:
      quantity() {
      }
      virtual ~quantity() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };
  class expression : public SAXObject, public IdentifiableExpressionType {
    public:
      expression() {
      }
      virtual ~expression() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };
  class position : public SAXObject, public positionType {
    public:
      position() {
      }
      virtual ~position() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };
  class rotation : public SAXObject, public rotationType {
    public:
      rotation() {
      }
      virtual ~rotation() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };
};

#endif // DEFINE_H
