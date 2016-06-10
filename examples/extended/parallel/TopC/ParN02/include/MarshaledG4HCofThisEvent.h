// This file was generated automatically by marshalgen.

#ifndef MarshaledG4HCofThisEvent_H
#define MarshaledG4HCofThisEvent_H


#include <G4HCofThisEvent.hh>
//MSH_include_begin
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "ExN02TrackerHit.hh"
#include "MarshaledExN02TrackerHit.h"
#include "MarshaledG4THitsCollection.h"
#include "MarshaledG4VHitsCollection.h"
//MSH_include_end

#include <stdio.h>
#include <string.h>
#include "MarshaledObj.h"

  class MarshaledG4HCofThisEvent;

  class ShadowedMarshaledG4HCofThisEvent : public G4HCofThisEvent{
    friend class MarshaledG4HCofThisEvent;
};

  class MarshaledG4HCofThisEvent : public MarshaledObj {
public:
    G4HCofThisEvent* param;
    ShadowedMarshaledG4HCofThisEvent* Shadowed_param;
public:


// Function implementations

MarshaledG4HCofThisEvent(G4HCofThisEvent* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledG4HCofThisEvent*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
}

MarshaledG4HCofThisEvent(void *buf, char chIsUnmarshaling = 'u')
: MarshaledObj(buf, chIsUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledG4HCofThisEvent() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

G4HCofThisEvent* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
        param = new G4HCofThisEvent();
        }
        this->Shadowed_param = (ShadowedMarshaledG4HCofThisEvent*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        return this->param;
    }
}

void unmarshalTo(G4HCofThisEvent* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledG4HCofThisEvent*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
    }
}

void marshal1() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); 
       // 4 bytes for the total size of field, 4 bytes for the number 
       // of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total 
    //size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
        int copy_off = 0;
        int elementNum;
         elementNum = param->GetNumberOfCollections(); 
        memcpy( msh_cursor+copy_off, &elementNum,sizeof(int));
        copy_off += sizeof(int);
        for(int index=0;index<elementNum;index++){
          G4VHitsCollection* anElement;
           anElement = param->GetHC(index); 
          MarshaledG4VHitsCollection marEle(anElement);
          EXTEND_BUFFER(marEle.getBufferSize());
          memcpy(msh_cursor+copy_off, marEle.getBuffer(), marEle.getBufferSize());
          copy_off += marEle.getBufferSize();
        }
        msh_currentSize = copy_off;

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
       int copy_off = 0;
       int elementNum;
       memcpy(&elementNum, msh_cursor+copy_off, sizeof(int));
       copy_off += sizeof(int);
       for(int index=0;index<elementNum;index++){
         MarshaledG4VHitsCollection marEle(msh_cursor+copy_off);
         G4VHitsCollection* anElement = (G4VHitsCollection*)marEle.unmarshal();
         copy_off += marEle.getBufferSize();
         param->AddHitsCollection(index, anElement); 
       }

    }
    msh_cursor += msh_currentSize;
}

};
#endif

