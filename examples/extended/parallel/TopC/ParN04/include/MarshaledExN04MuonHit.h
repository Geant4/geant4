// This file was generated automatically by marshalgen.

#ifndef MarshaledExN04MuonHit_H
#define MarshaledExN04MuonHit_H


#include <ExN04MuonHit.hh>
//MSH_include_begin
#include "MarshaledG4String.h"
//MSH_include_end

#include <stdio.h>
#include <string.h>
#include "MarshaledObj.h"

  class MarshaledExN04MuonHit;

  class ShadowedMarshaledExN04MuonHit : public ExN04MuonHit{
    friend class MarshaledExN04MuonHit;
};

  class MarshaledExN04MuonHit : public MarshaledObj {
public:
    ExN04MuonHit* param;
    ShadowedMarshaledExN04MuonHit* Shadowed_param;
public:


// Function implementations

MarshaledExN04MuonHit(ExN04MuonHit* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledExN04MuonHit*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
    marshal2();
}

MarshaledExN04MuonHit(void *buf, char chIsUnmarshaling = 'u')
: MarshaledObj(buf, chIsUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledExN04MuonHit() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

ExN04MuonHit* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
        param = new ExN04MuonHit();
        }
        this->Shadowed_param = (ShadowedMarshaledExN04MuonHit*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
        return this->param;
    }
}

void unmarshalTo(ExN04MuonHit* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledExN04MuonHit*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
    }
}

void marshal1() {
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
    //of this field here later
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

void unmarshal1() {
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

void marshal2() {
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

void unmarshal2() {
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

