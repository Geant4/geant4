#ifndef IDENTIFIABLEVOLUMETYPE_H
#define IDENTIFIABLEVOLUMETYPE_H 1

#include <string>

class IdentifiableVolumeType {
public:
  IdentifiableVolumeType() {
  }
  ~IdentifiableVolumeType() {
  }
  const std::string& get_name() const {
    return m_ID;
  }
  void set_name( const std::string& n ) {
    m_ID = n;
  }
private:
  std::string m_ID;
};



#endif // IDENTIFIABLEVOLUMETYPE_H
