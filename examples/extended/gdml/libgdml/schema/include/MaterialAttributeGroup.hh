#ifndef MATERIALATTRIBUTEGROUP_H
#define MATERIALATTRIBUTEGROUP_H 1

#include <string>

/*
enum MaterialState {
  unknown = 0,
  gas,
  solid,
  liquid
};
*/

struct MaterialAttributeGroup {
  std::string m_ID;
  std::string m_formula;
  //MaterialState m_state;
  std::string m_state;
};

#endif // MATERIALATTRIBUTEGROUP_H
