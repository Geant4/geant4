#ifndef EXPRESSIONTYPE_H
#define EXPRESSIONTYPE_H 1


class ExpressionType {
public:
  ExpressionType() {
  }
  ~ExpressionType() {
  }
  void set_text( std::string& s ) {
    m_text = s;
  }
  void set_text( char* s ) {
    m_text = s;
  }
  std::string get_text() const {
    return m_text;
  }
private:
  // mixed
  std::string m_text;
};


#endif // EXPRESSIONTYPE_H
