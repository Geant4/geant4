#ifndef G4FUNCTORIDENTIFIER_HH
#define G4FUNCTORIDENTIFIER_HH

/*
struct G4FunctorIdentifier : public G4String {
  static G4String GetDelimieter() {return ":";}
};
*/

class G4FunctorIdentifier {

public:
  typedef unsigned long Key;

  G4FunctorIdentifier(const G4String& name = "null",
		      Key parentID=0)
  {
    fParentID = parentID;
    fUniqueID = NextKey();
    fName = name;
  }

  Key UniqueID() const {return fUniqueID;}
  Key ParentID() const {return fParentID;}
  const G4String& Name() const {return fName;}

  bool operator == (const G4FunctorIdentifier& other) const
  {
    return ((fName == other.fName) &&
	    (fUniqueID == other.fUniqueID) &&
	    (fParentID == other.fParentID));
  }
private:

  Key NextKey() const {
    static Key nKey = 0;
    return ++nKey;
  }

  G4String fName;
  Key fUniqueID;
  Key fParentID;

};

#endif
