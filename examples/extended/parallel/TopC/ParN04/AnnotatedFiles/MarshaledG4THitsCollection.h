//
/// \file MarshaledG4THitsCollection.h
/// \brief Definition of the MaraledG4THitsCollection class
//
// This file was generated automatically by marshalgen.

#ifndef MarshaledG4HitsCollection_H
#define MarshaledG4HitsCollection_H


#include <G4THitsCollection.hh>
//MSH_include_begin
#include "ExN04CalorimeterHit.hh"
#include "ExN04MuonHit.hh"
#include "ExN04TrackerHit.hh"
#include "MarshaledExN04CalorimeterHit.h"
#include "MarshaledExN04MuonHit.h"
#include "MarshaledExN04TrackerHit.h"
//MSH_include_end

#include <stdio.h>
#include <string.h>
#include "MarshaledObj.h"

  class MarshaledG4HitsCollection;

  class ShadowedMarshaledG4HitsCollection : public G4HitsCollection{
    friend class MarshaledG4HitsCollection;
};


template <class T> class MarshaledG4THitsCollection;

template <class T> class ShadowedMarshaledG4THitsCollection : public G4THitsCollection<T>{
    friend class MarshaledG4THitsCollection<T>;
};

  class MarshaledG4HitsCollection : public MarshaledObj {
public:
    G4HitsCollection* param;
    ShadowedMarshaledG4HitsCollection* Shadowed_param;
public:


// Function implementations

MarshaledG4HitsCollection(G4HitsCollection* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledG4HitsCollection*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
}

MarshaledG4HitsCollection(void *buf, char isUnmarshaling = 'u')
: MarshaledObj(buf, isUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledG4HitsCollection() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

G4HitsCollection* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
	param = new G4HitsCollection();
	}
        this->Shadowed_param = (ShadowedMarshaledG4HitsCollection*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        return this->param;
    }
}

void unmarshalTo(G4HitsCollection* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledG4HitsCollection*)this->param;
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
    EXTEND_BUFFER(msh_currentSize + sizeof(int) + sizeof(int)); // 4 bytes for the total size of field, 4 bytes for the number of elements in the array (in the case of array marshaling)
    //Mark the beginning position for this field, will write the total size of this field here later
    msh_field_begin = msh_cursor;

    //Advance cursor of distance = sizeof(int)
    msh_cursor += sizeof(int);

    //Now just copy "get" functions here
    {
	int copy_off = 0;
	int elementNum;
	 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        elementNum = ((G4THitsCollection<ExN04CalorimeterHit>*)param)->entries(); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        elementNum = ((G4THitsCollection<ExN04MuonHit>*)param)->entries(); 
     else
        elementNum = ((G4THitsCollection<ExN04TrackerHit>*)param)->entries(); 
    
	memcpy( msh_cursor+copy_off, &elementNum,sizeof(int));
	copy_off += sizeof(int);
	for(int index=0;index<elementNum;index++){
		if(0){}
		else if((dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL) ){
			void* anElement;
			 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04CalorimeterHit>*)param))[index]; 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04MuonHit>*)param))[index]; 
     else
        anElement = (*((G4THitsCollection<ExN04TrackerHit>*)param))[index]; 
    
			MarshaledExN04CalorimeterHit  marEle((ExN04CalorimeterHit* )anElement);
			EXTEND_BUFFER(marEle.getBufferSize());
			memcpy(msh_cursor+copy_off, marEle.getBuffer(), marEle.getBufferSize());
			copy_off += marEle.getBufferSize();
		}
		else if( (dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL) ){
			void* anElement;
			 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04CalorimeterHit>*)param))[index]; 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04MuonHit>*)param))[index]; 
     else
        anElement = (*((G4THitsCollection<ExN04TrackerHit>*)param))[index]; 
    
			MarshaledExN04MuonHit  marEle((ExN04MuonHit* )anElement);
			EXTEND_BUFFER(marEle.getBufferSize());
			memcpy(msh_cursor+copy_off, marEle.getBuffer(), marEle.getBufferSize());
			copy_off += marEle.getBufferSize();
		}
		else if( true ){
			void* anElement;
			 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04CalorimeterHit>*)param))[index]; 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        anElement = (*((G4THitsCollection<ExN04MuonHit>*)param))[index]; 
     else
        anElement = (*((G4THitsCollection<ExN04TrackerHit>*)param))[index]; 
    
			MarshaledExN04TrackerHit marEle((ExN04TrackerHit*)anElement);
			EXTEND_BUFFER(marEle.getBufferSize());
			memcpy(msh_cursor+copy_off, marEle.getBuffer(), marEle.getBufferSize());
			copy_off += marEle.getBufferSize();
		}
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
	if(0){}
	else if((dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL) ){
		int copy_off = 0;
		int elementNum;
		memcpy(&elementNum, msh_cursor+copy_off, sizeof(int));
		copy_off += sizeof(int);
		for(int index=0;index<elementNum;index++){
			MarshaledExN04CalorimeterHit  marEle(msh_cursor+copy_off);
			ExN04CalorimeterHit*  anElement = (ExN04CalorimeterHit* )marEle.unmarshal();
			copy_off += marEle.getBufferSize();
			
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04CalorimeterHit>*)param)->insert((ExN04CalorimeterHit*)anElement); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04MuonHit>*)param)->insert((ExN04MuonHit*)anElement); 
     else
        ((G4THitsCollection<ExN04TrackerHit>*)param)->insert((ExN04TrackerHit*)anElement); 
    
		}
	}
	else if( (dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL) ){
		int copy_off = 0;
		int elementNum;
		memcpy(&elementNum, msh_cursor+copy_off, sizeof(int));
		copy_off += sizeof(int);
		for(int index=0;index<elementNum;index++){
			MarshaledExN04MuonHit  marEle(msh_cursor+copy_off);
			ExN04MuonHit*  anElement = (ExN04MuonHit* )marEle.unmarshal();
			copy_off += marEle.getBufferSize();
			
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04CalorimeterHit>*)param)->insert((ExN04CalorimeterHit*)anElement); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04MuonHit>*)param)->insert((ExN04MuonHit*)anElement); 
     else
        ((G4THitsCollection<ExN04TrackerHit>*)param)->insert((ExN04TrackerHit*)anElement); 
    
		}
	}
	else if( true ){
		int copy_off = 0;
		int elementNum;
		memcpy(&elementNum, msh_cursor+copy_off, sizeof(int));
		copy_off += sizeof(int);
		for(int index=0;index<elementNum;index++){
			MarshaledExN04TrackerHit marEle(msh_cursor+copy_off);
			ExN04TrackerHit* anElement = (ExN04TrackerHit*)marEle.unmarshal();
			copy_off += marEle.getBufferSize();
			
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04CalorimeterHit>*)param)->insert((ExN04CalorimeterHit*)anElement); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>(param)!=NULL)
        ((G4THitsCollection<ExN04MuonHit>*)param)->insert((ExN04MuonHit*)anElement); 
     else
        ((G4THitsCollection<ExN04TrackerHit>*)param)->insert((ExN04TrackerHit*)anElement); 
    
		}
	}

    }
    msh_cursor += msh_currentSize;
}

};
template <class T> class MarshaledG4THitsCollection : public MarshaledObj {
public:
    G4THitsCollection<T>* param;
    ShadowedMarshaledG4THitsCollection<T>* Shadowed_param;
public:


// Function implementations

MarshaledG4THitsCollection(G4THitsCollection<T>* objptr) : MarshaledObj() {
    msh_isUnmarshalDone = false;
    this->param = objptr;
    this->Shadowed_param = (ShadowedMarshaledG4THitsCollection<T>*)this->param;
    if (objptr == NULL)
        return;

    marshal1();
}

MarshaledG4THitsCollection(void *buf, char isUnmarshaling = 'u')
: MarshaledObj(buf, isUnmarshaling) {
    msh_isUnmarshalDone = false;
}

~MarshaledG4THitsCollection() {
    //if(msh_isUnmarshalDone && this->param != NULL) {
        //delete this->param;
    //}
}

G4THitsCollection<T>* unmarshal() {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return NULL;
    } else {
        {
	param = new G4THitsCollection<T>();
	}
        this->Shadowed_param = (ShadowedMarshaledG4THitsCollection<T>*)this->param;
        this->msh_isUnmarshalDone = true;
        unmarshal1();
        return this->param;
    }
}

void unmarshalTo(G4THitsCollection<T>* obj) {
    //We don't want to unmarshal the buffer is empty.
    if(msh_size <= MSH_HEADER_SIZE) {
        //This is buggy, we can't always assume that
        //obj == NULL <==> List is empty.
        return;
    } else {
        this->param = obj;
        this->Shadowed_param = (ShadowedMarshaledG4THitsCollection<T>*)this->param;
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
		MarshaledG4HitsCollection marParent(param);
		EXTEND_BUFFER(marParent.getBufferSize());
		memcpy(msh_cursor,marParent.getBuffer(), marParent.getBufferSize());
		msh_currentSize = marParent.getBufferSize();

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
		MarshaledG4HitsCollection marObj(msh_cursor);
		marObj.unmarshalTo(param);

    }
    msh_cursor += msh_currentSize;
}

};
#endif

