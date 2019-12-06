/// \file MarshaledG4VHitsCollection.h
/// \brief Definition of the MaraledG4VHitsCollection class
//
// This file was generated automatically by marshalgen.

#ifndef MarshaledG4VHitsCollection_H
#define MarshaledG4VHitsCollection_H


#include <G4VHitsCollection.hh>


#include <stdio.h>
#include <string.h>
#include "MarshaledObj.h"

  class MarshaledG4VHitsCollection;

  class ShadowedMarshaledG4VHitsCollection : public G4VHitsCollection{
    friend class MarshaledG4VHitsCollection;
};

  class MarshaledG4VHitsCollection : public MarshaledObj {
public:
    G4VHitsCollection* param;
    ShadowedMarshaledG4VHitsCollection* Shadowed_param;
public:


// Function implementations

MarshaledG4VHitsCollection(G4VHitsCollection* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledG4VHitsCollection*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
    marshal2();
    marshal3();
}

MarshaledG4VHitsCollection(void *buf, char isUnmarshaling = 'u')
: MarshaledObj(buf, isUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledG4VHitsCollection() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

G4VHitsCollection* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
	if(0){}
	else if(msh_typechoice == 0){
	param = new G4THitsCollection<ExN04CalorimeterHit>("","");
	}
	else if(msh_typechoice == 1){
	param = new G4THitsCollection<ExN04MuonHit>("","");
	}
	else if(msh_typechoice == 2){
	param = new G4THitsCollection<ExN04TrackerHit>("","");
	}
	}
        this->Shadowed_param = (ShadowedMarshaledG4VHitsCollection*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
        unmarshal3();
        return this->param;
    }
}

void unmarshalTo(G4VHitsCollection* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledG4VHitsCollection*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        unmarshal2();
        unmarshal3();
    }
}

void marshal1() {
    //declare field_size to be the size of this field
    int msh_currentSize = 0;
    if (isUnmarshaling())
        throw "Tried to marshal in obj marked isUnmarshaling == true";

    //Copy the sizespec into msh_currentSize here:
    {
	// no need to declare size since msh_currentSize is already assigned in the MARSHAL field

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); // 4 bytes for the total size of field, 4 bytes for the number of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
	G4String anElement;
	 anElement = param->GetName(); 
	MarshaledG4String var(&anElement);
	EXTEND_BUFFER(var.getBufferSize());
	msh_currentSize = var.getBufferSize();
	memcpy(msh_cursor, var.getBuffer(), var.getBufferSize());
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
	MarshaledG4String var(msh_cursor, 'u');
	G4String anElement;
	var.unmarshalTo(&anElement);
	 Shadowed_param->collectionName=anElement; 

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
	// no need to declare size since msh_currentSize is already assigned in the MARSHAL field

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); // 4 bytes for the total size of field, 4 bytes for the number of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
	G4String anElement;
	 anElement = param->GetSDname(); 
	MarshaledG4String var(&anElement);
	EXTEND_BUFFER(var.getBufferSize());
	msh_currentSize = var.getBufferSize();
	memcpy(msh_cursor, var.getBuffer(), var.getBufferSize());
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
	MarshaledG4String var(msh_cursor, 'u');
	G4String anElement;
	var.unmarshalTo(&anElement);
	 Shadowed_param->SDname=anElement; 

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
		//code for size, just dummy code because the size will be set correctly at the end of marshaling code

    }

    //Increase the size of buffer if needed
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); // 4 bytes for the total size of field, 4 bytes for the number of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
	if(0){}
	else if((param->GetName() == "calCollection") ){
		G4THitsCollection<ExN04CalorimeterHit> *aObj621 = (G4THitsCollection<ExN04CalorimeterHit>*)param;
		MarshaledG4THitsCollection<ExN04CalorimeterHit> marChild(aObj621);
		EXTEND_BUFFER(marChild.getBufferSize());
		memcpy(msh_cursor,marChild.getBuffer(), marChild.getBufferSize());
		msh_currentSize = marChild.getBufferSize();
		msh_typechoice = 0;
	}
	else if( (param->GetName() == "muonCollection") ){
		G4THitsCollection<ExN04MuonHit> *aObj621 = (G4THitsCollection<ExN04MuonHit>*)param;
		MarshaledG4THitsCollection<ExN04MuonHit> marChild(aObj621);
		EXTEND_BUFFER(marChild.getBufferSize());
		memcpy(msh_cursor,marChild.getBuffer(), marChild.getBufferSize());
		msh_currentSize = marChild.getBufferSize();
		msh_typechoice = 1;
	}
	else if( true ){
		G4THitsCollection<ExN04TrackerHit> *aObj621 = (G4THitsCollection<ExN04TrackerHit>*)param;
		MarshaledG4THitsCollection<ExN04TrackerHit> marChild(aObj621);
		EXTEND_BUFFER(marChild.getBufferSize());
		memcpy(msh_cursor,marChild.getBuffer(), marChild.getBufferSize());
		msh_currentSize = marChild.getBufferSize();
		msh_typechoice = 2;
	}

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
	if(0){}
	else if(msh_typechoice == 0){
		MarshaledG4THitsCollection<ExN04CalorimeterHit> marObj(msh_cursor);
		marObj.unmarshalTo((G4THitsCollection<ExN04CalorimeterHit>*)param);
	}
	else if(msh_typechoice == 1){
		MarshaledG4THitsCollection<ExN04MuonHit> marObj(msh_cursor);
		marObj.unmarshalTo((G4THitsCollection<ExN04MuonHit>*)param);
	}
	else if(msh_typechoice == 2){
		MarshaledG4THitsCollection<ExN04TrackerHit> marObj(msh_cursor);
		marObj.unmarshalTo((G4THitsCollection<ExN04TrackerHit>*)param);
	}

    }
    msh_cursor += msh_currentSize;
}

};
#endif

