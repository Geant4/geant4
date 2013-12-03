// This file was generated automatically by marshalgen.

#ifndef MarshaledExN04CalorimeterHit_H
#define MarshaledExN04CalorimeterHit_H


#include <ExN04CalorimeterHit.hh>
//MSH_include_begin
#include "MarshaledG4String.h"
//MSH_include_end

#include <stdio.h>
#include <string.h>
#include "MarshaledObj.h"

  class MarshaledExN04CalorimeterHit;

  class ShadowedMarshaledExN04CalorimeterHit : public ExN04CalorimeterHit{
    friend class MarshaledExN04CalorimeterHit;
};

  class MarshaledExN04CalorimeterHit : public MarshaledObj {
public:
    ExN04CalorimeterHit* param;
    ShadowedMarshaledExN04CalorimeterHit* Shadowed_param;
public:


// Function implementations

MarshaledExN04CalorimeterHit(ExN04CalorimeterHit* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledExN04CalorimeterHit*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
    marshal2();
    marshal3();
    marshal4();
}

MarshaledExN04CalorimeterHit(void *buf, char chIsUnmarshaling = 'u')
: MarshaledObj(buf, chIsUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledExN04CalorimeterHit() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

ExN04CalorimeterHit* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
        param = new ExN04CalorimeterHit();
        }
        this->Shadowed_param = (ShadowedMarshaledExN04CalorimeterHit*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
        unmarshal3();
        unmarshal4();
        return this->param;
    }
}

void unmarshalTo(ExN04CalorimeterHit* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledExN04CalorimeterHit*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
        unmarshal3();
        unmarshal4();
    }
}

void marshal1() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {
        msh_currentSize = sizeof(G4int);

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); 
       // 4 bytes for the total size of field, 4 bytes for the number 
       // of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total
    // size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
        G4int anElement;
         anElement = param->GetZ(); 
        memcpy(msh_cursor, &anElement, sizeof(G4int));
    }
    //Now advance the cursor
    msh_cursor += msh_currentSize;
    //Now set the size of this field
    int tmp; //use memcpy instead of *(int*)... =... to prevent bus error
    tmp = (msh_cursor-msh_field_begin) - sizeof(int);
    memcpy(msh_field_begin, &tmp, sizeof(int));

    //Now set msh_size
    msh_size = msh_cursor - msh_buffer;
    MSH_SET_TOTALSIZE(msh_size);    MSH_SET_TYPECHOICE(msh_typechoice);
}

void unmarshal1() {
    //declare currentSize to be the size of this field
    int msh_currentSize = 0;
    //copy the size of the current field into currentSize
    memcpy(&msh_currentSize, msh_cursor, sizeof(int));
    msh_cursor += sizeof(int);
    //Now copy the setspec here
    {
        G4int anElement;
        memcpy(&anElement, msh_cursor, sizeof(G4int));
         param->SetCellID(anElement, param->GetPhi()); 

    }
    msh_cursor += msh_currentSize;
}

void marshal2() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {
        msh_currentSize = sizeof(G4int);

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); 
       // 4 bytes for the total size of field, 4 bytes for the number of 
       // elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total 
    //size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
        G4int anElement;
         anElement = param->GetPhi(); 
        memcpy(msh_cursor, &anElement, sizeof(G4int));
    }
    //Now advance the cursor
    msh_cursor += msh_currentSize;
    //Now set the size of this field
    int tmp; //use memcpy instead of *(int*)... =... to prevent bus error
    tmp = (msh_cursor-msh_field_begin) - sizeof(int);
    memcpy(msh_field_begin, &tmp, sizeof(int));

    //Now set msh_size
    msh_size = msh_cursor - msh_buffer;
    MSH_SET_TOTALSIZE(msh_size);    MSH_SET_TYPECHOICE(msh_typechoice);
}

void unmarshal2() {
    //declare currentSize to be the size of this field
    int msh_currentSize = 0;
    //copy the size of the current field into currentSize
    memcpy(&msh_currentSize, msh_cursor, sizeof(int));
    msh_cursor += sizeof(int);
    //Now copy the setspec here
    {
        G4int anElement;
        memcpy(&anElement, msh_cursor, sizeof(G4int));
         param->SetCellID(param->GetZ(), anElement); 

    }
    msh_cursor += msh_currentSize;
}

void marshal3() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {
        msh_currentSize = sizeof(G4double);

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); 
       // 4 bytes for the total size of field, 4 bytes for the number of 
       // elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size
    // of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
        G4double anElement;
         anElement = param->GetEdep(); 
        memcpy(msh_cursor, &anElement, sizeof(G4double));
    }
    //Now advance the cursor
    msh_cursor += msh_currentSize;
    //Now set the size of this field
    int tmp; //use memcpy instead of *(int*)... =... to prevent bus error
    tmp = (msh_cursor-msh_field_begin) - sizeof(int);
    memcpy(msh_field_begin, &tmp, sizeof(int));

    //Now set msh_size
    msh_size = msh_cursor - msh_buffer;
    MSH_SET_TOTALSIZE(msh_size);    MSH_SET_TYPECHOICE(msh_typechoice);
}

void unmarshal3() {
    //declare currentSize to be the size of this field
    int msh_currentSize = 0;
    //copy the size of the current field into currentSize
    memcpy(&msh_currentSize, msh_cursor, sizeof(int));
    msh_cursor += sizeof(int);
    //Now copy the setspec here
    {
        G4double anElement;
        memcpy(&anElement, msh_cursor, sizeof(G4double));
         param->SetEdep(anElement); 

    }
    msh_cursor += msh_currentSize;
}

void marshal4() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {
        msh_currentSize = sizeof(G4ThreeVector);

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); 
       // 4 bytes for the total size of field, 4 bytes for the number of
       // elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size
    //of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
        G4ThreeVector anElement;
         anElement = param->GetPos(); 
        memcpy(msh_cursor, &anElement, sizeof(G4ThreeVector));
    }
    //Now advance the cursor
    msh_cursor += msh_currentSize;
    //Now set the size of this field
    int tmp; //use memcpy instead of *(int*)... =... to prevent bus error
    tmp = (msh_cursor-msh_field_begin) - sizeof(int);
    memcpy(msh_field_begin, &tmp, sizeof(int));

    //Now set msh_size
    msh_size = msh_cursor - msh_buffer;
    MSH_SET_TOTALSIZE(msh_size);    MSH_SET_TYPECHOICE(msh_typechoice);
}

void unmarshal4() {
    //declare currentSize to be the size of this field
    int msh_currentSize = 0;
    //copy the size of the current field into currentSize
    memcpy(&msh_currentSize, msh_cursor, sizeof(int));
    msh_cursor += sizeof(int);
    //Now copy the setspec here
    {
        G4ThreeVector anElement;
        memcpy(&anElement, msh_cursor, sizeof(G4ThreeVector));
         param->SetPos(anElement); 

    }
    msh_cursor += msh_currentSize;
}

};
#endif

