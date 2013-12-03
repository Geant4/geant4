 /**********************************************************************
  *                Include file with Base Class of Marshalgen          *
  **********************************************************************/

#ifndef MARSHALEDOBJ_H
#define MARSHALEDOBJ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MSH_ASSERT(X) {assert(X);}

#define MSH_HEADER_SIZE (sizeof(int)*2)
// the first field (first sizeof(int) bytes) contains the $TYPE_CHOICE
// the second field contains the size including the header
#define MSH_TOTALSIZE_OFFSET (sizeof(int))
#define MSH_TYPECHOICE_OFFSET 0

#define MSH_SET_TYPECHOICE(X) { memcpy(msh_buffer+MSH_TYPECHOICE_OFFSET,&(X),sizeof(int));}
#define MSH_SET_TOTALSIZE(X) { memcpy(msh_buffer+MSH_TOTALSIZE_OFFSET,&(X),sizeof(int));}
#define MSH_GET_TYPECHOICE(X,BUF) { memcpy(&(X), ((char*)BUF)+MSH_TYPECHOICE_OFFSET,sizeof(int));}
#define MSH_GET_TOTALSIZE(X,BUF) { memcpy(&(X), ((char*)BUF)+MSH_TOTALSIZE_OFFSET,sizeof(int));}


class MarshaledObj {
  private:
    // Make sure all marshaled objects are word aligned.
    static const int WORD_SIZE = sizeof(long);
  public:
    static  int ROUND_UP( int x ){
                return (((x)+(WORD_SIZE-1)) / WORD_SIZE) * WORD_SIZE;
        }

  public:
        // Constructs an empty MarshaledObj, 
        MarshaledObj(){
                msh_extent = 128;
                msh_size = MSH_HEADER_SIZE;
                msh_isUnmarshalDone = false;

                msh_buffer = (char *)malloc(msh_extent);
                MSH_ASSERT(msh_buffer);

                msh_cursor = msh_buffer + MSH_HEADER_SIZE;
                msh_field_begin = msh_cursor;

                msh_typechoice = 0;
                int totalsize = msh_cursor-msh_buffer;

                MSH_SET_TYPECHOICE(msh_typechoice);
                MSH_SET_TOTALSIZE(totalsize);
        }

    //MarshaledObj(void *buf);
        // This constructs a MarshledObj from a buffer (of type char*) for unmarshaling.
        // buf is obtain from an already marshaled object.
        // The first field of buf must be an int that contains the size of the buf
        // NOT including itself.
        // isUnmarshaling must be 'u' (for unmarshaling) .
    MarshaledObj(void *buf, char chIsUnmarshaling) {
                msh_isUnmarshalDone = false;

                if(chIsUnmarshaling != 'u') {
                        printf("MarshaledObj(void*, char): wrong argument\n");
                        return;
                }

                //msh_extent = ROUND_UP(*(int *)buf + sizeof(int));
                MSH_GET_TYPECHOICE(msh_typechoice,buf);
                                        
                MSH_GET_TOTALSIZE(msh_size,buf);
                msh_extent = ROUND_UP(msh_size);

                msh_buffer = (char *)malloc(msh_extent);
                MSH_ASSERT(msh_buffer);

                memcpy(msh_buffer, (char *)buf, msh_extent);
                msh_cursor = msh_buffer + MSH_HEADER_SIZE;
                msh_field_begin = msh_cursor;

                //MSH_SET_TYPECHOICE(msh_typechoice);

        }

    ~MarshaledObj() {
                if ( ! isUnmarshaling() )
                        free(msh_buffer);
        }

    inline  bool isUnmarshaling() {
                return (msh_extent <= 0);
        }

  private:
    // Dont use copy constructor
    const MarshaledObj& operator=(const MarshaledObj& right);

  protected:
    int msh_typechoice;  // alias of $TYPE_CHOICE

    // points to the buffer (header+body)                                
    char *msh_buffer;

    // msh_field_begin points to the size of the current field being marshaled
    char* msh_field_begin;

    // msh_size contains the total size of msh_buffer. i.e.,
    size_t msh_size;

    // msh_cursor points to the next field to be marshaled.
    char *msh_cursor;

    // msh_extent is the total allocated space for msh_buffer.
    // msh_extent is always >= msh_size
    size_t msh_extent;

    bool msh_isUnmarshalDone; //Is unmarshaling done yet?

  public:
    inline void EXTEND_BUFFER(int size){
                msh_size += size;
                if(msh_size > msh_extent){
                        resizeBuffer(msh_size);
                }
    }

        void resizeBuffer(size_t new_size ) {
                int msh_displacement = msh_cursor - msh_buffer;
                int field_displacement = msh_field_begin - msh_buffer;

                while(new_size > msh_extent)
                        msh_extent *= 2;

                msh_buffer = (char *)realloc( msh_buffer, msh_extent);
                MSH_ASSERT(msh_buffer);

                msh_cursor = msh_buffer + msh_displacement;
                msh_field_begin = msh_buffer + field_displacement;
        }

  public:
    // Returns the total size of buffer
    inline int getBufferSize() {
        return msh_size;
    }

    inline char *getBuffer() {
        return msh_buffer;
    }

    /* p: pointer to the data field, size: size of that primitive data field */
    void marshalPrimitive(void* p, int size) {
                int msh_currentSize;
                if (isUnmarshaling())
                        throw "Tried to marshal in object marked isUnmarshaling = true";
                msh_currentSize =  size;
                EXTEND_BUFFER(msh_currentSize + sizeof(int));

                // *(int *)msh_cursor = msh_currentSize;
                memcpy(msh_cursor, &msh_currentSize, sizeof(int));
                msh_cursor += sizeof(int);
                memcpy(msh_cursor, p, size);
                msh_cursor += msh_currentSize;
                msh_size = msh_cursor - msh_buffer;

                MSH_SET_TOTALSIZE(msh_size);
    }

    void unmarshalPrimitive(void* p, int size) {
                int msh_currentSize;
                //memcpy(&msh_currentSize, msh_cursor, sizeof(int));
                /* in case *msh_cursor is invalid, use "size" not to crash the memory */
                msh_currentSize = size;
                msh_cursor += sizeof(int);
                memcpy(p, msh_cursor, msh_currentSize);
                msh_cursor += msh_currentSize;
                //msh_size = msh_cursor - msh_buffer;
    }
};

/* Used for distinguish the class types of the template parameter 
   vietha 2003.05.01 */
template <class T,class>
class MSH_IsSameClass
{
public:
  enum {Is = 0};
};

template<class T>
class MSH_IsSameClass<T,T>
{
public:
  enum {Is = 1};
};

#endif
