#ifndef _InputReader_H
#define _InputReader_H

#ifdef IS_GCC
#pragma interface
#endif
              
#include "DList.hh"
#include "String.hh"
#include "Error.hh"

class InputItem;

class InputReader {
  public:                
    InputReader();
    ~InputReader() {} 
    void setCommentLeader(const String& cl);
    InputReader& operator>>(InputItem& item);
    friend istream& operator>>(istream& is, InputReader& reader);
    
    class ErrDuplicateItem 
      : public Error 
    {
      public:
        ErrDuplicateItem(const String& key);
        ErrDuplicateItem(const ErrDuplicateItem& );
        ~ErrDuplicateItem() {}
        virtual void writeMessage(ostream& os) const;
      private:      
        String theKey;
    };
    
    class ErrDuplicateInput 
      : public Error 
    {
      public:
        ErrDuplicateInput(const String& key, int count);
        ErrDuplicateInput(const ErrDuplicateInput&);
        ~ErrDuplicateInput() {}
        virtual void writeMessage(ostream& os) const;
      private:      
        String theKey;
	int theCount;
    };
    
    class ErrReadFailed
      : public Error 
    {
      public:
        ErrReadFailed(const String& key, int count);
        ErrReadFailed(const ErrReadFailed& err);
        ~ErrReadFailed() {}
        virtual void writeMessage(ostream& os) const;
      private:      
        String theKey;
	int theCount;
    };

    class ErrNotKnown
      : public Error 
    {
      public:
        ErrNotKnown(const String& key, int count);
        ErrNotKnown(const ErrNotKnown& err);
        ~ErrNotKnown() {}
        virtual void writeMessage(ostream& os) const;
      private:      
        String theKey;
	int theCount;
    };
  private:
    DList<InputItem> itemList;
    String commentLeader;
};

#endif // _InputReader_H

