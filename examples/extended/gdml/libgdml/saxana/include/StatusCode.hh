#ifndef STATUS_CODE_H
#define STATUS_CODE_H 1

class StatusCode{
public:
  typedef enum
  {
    eError = 0,
    eOk = 1,
    eWarning,
    eParserInitFailure
  } EInfo;

public:
  StatusCode( EInfo info = StatusCode::eOk )
  : fInfo( info )
  {
  }  
  ~StatusCode()
  {
  }
  
  bool IsOk()
  {
    return( (fInfo == eOk) ? true : false );
  } 
  bool IsFailure()
  {
    return( (fInfo != eOk) ? true : false );
  }
  StatusCode::EInfo Info()
  {
    return fInfo;
  }
  
private:
  EInfo fInfo;
};

#endif // STATUS_CODE_H

