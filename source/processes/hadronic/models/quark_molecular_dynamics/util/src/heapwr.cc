#include "heapwr.hh"

//
// Funktionen, die den Heap-Wrapper implementieren
//
//
// werror - Drucke einen Heap-Fehler aus und beende Programm
//

Wrapper *pFirst = 0;  // Pointer auf die verkettete Liste der
                      // allokierten Heap-Blocks

void werror (Wrapper *pBlock, char *pErrorString, int fatal) {
    printf("Heap error: %s, Allocated: %s, %d\n", pErrorString,
            pBlock->fileName, pBlock->lineNumber);
    if (fatal) exit(1);
}

//
// new - allokiere einen Block vom Heap und setze einen tag auf den
//       Debugging-Typ und die Dateninformation - gib einen Pointer
//       auf die User-Daten zur"uck
//
#ifndef __CLEAR__
void* operator new(size_t sizeBlock, char *pFile, int lineNo) {
    return allocateFn(sizeBlock, pFile, lineNo);
}

void* operator new(size_t sizeBlock) {
    return allocateFn(sizeBlock, "Unknown", 0);
}
#ifdef IS_GCC
void* operator new[](size_t sizeBlock) {
    return allocateFn(sizeBlock, "Unknown", 0);
}
void* operator new[](size_t sizeBlock, char *pFile, int lineNo) {
    return allocateFn(sizeBlock, pFile, lineNo);
}
#endif
#endif
void* allocateFn(size_t sizeBlock, char *pFile, int lineNo) {
    Wrapper* pBlock;
    char* pUserData;
    if (sizeBlock == 0) {
        return (void*)0;
    } else {
        // nimm einen Block vom Heap
        pBlock = (Wrapper*) malloc (sizeBlock + sizeof(Wrapper));
        if (!pBlock)  {
            return (void*)0;
        }

	// speichere nun die Daten
	int len = strlen(pFile);
        int offset = (len >= NAMELENGTH) ? len-NAMELENGTH+1 : 0;
	strncpy(pBlock->fileName, pFile+offset, NAMELENGTH-1);
        pBlock->fileName[NAMELENGTH-1] = '\0';
        pBlock->lineNumber = lineNo;
        pBlock->length = sizeBlock;

        // speichere nun Epi- und Prolog
        pBlock->prologue = 0x12;
        pUserData = pBlock->data;
        char *pEpilogue = pUserData + sizeBlock;
        *pEpilogue = 0x21;

        // h"ange es an die Liste an
        pBlock->pNext = pFirst;
        pFirst = pBlock;
        // und gib dem User einen Zeiger auf seine Daten zur"uck
        return (void*)pUserData;
    }
}

//
// delete - akzeptiere jede Form der delete-Funktion (size-Information
//          ist unn"otig
//

#ifndef __CLEAR__
void operator delete(void* pUserData)  {
    deleteFn(pUserData);
}
#endif

//
// deleteFn - teste den Heap-Block um seine Validit"at zu pr"ufen,
//            bevor er aus der Liste gestrichen wird
//
void deleteFn(void* pUserData) {
#ifndef __CLEAR__
    // Berechne die Adresse des originalen Wrappers
    int offset = (int)&(((Wrapper*)0)->data);
    // (offset ist gerade die Anzahl von Bytes, die das Datenfeld
    //  vom Anfang des Wrapperblocks entfernt ist)
    Wrapper *pBlock = (Wrapper*)((char*)pUserData - offset);

    // Teste Prolog und Epilog
    if (pBlock->prologue != 0x12)
        werror (pBlock, "Prologue overwritten", 1);
    char *pEpilogue = (char*) pUserData + pBlock->length;
    if (*pEpilogue != 0x21)
        werror (pBlock, "Epilogue overwritten", 1);

    // Entferne den Block aus der Liste
    Wrapper *pWrapper = pFirst;
    int foundIt = 0;
    if (pWrapper == pBlock)  {
        pFirst = pBlock->pNext;
        foundIt = 1;
    }
    else {
        while (pWrapper)  {
            if (pWrapper->pNext == pBlock)  {
                pWrapper->pNext = pBlock->pNext;
                foundIt = 1;
                break;
            }
            pWrapper = pWrapper->pNext;
        }
    }
    if (!foundIt)  {
        werror (pBlock, "Block not in list", 1);
    }

    // nun gib den Block frei
    if (foundIt) free((void*)pBlock);
#else
    delete pUserData;
#endif
}

//
// displayAllocated - diese Funktion zeigt alle Blocks, die im Heap
//                    noch allokiert sind
//
void displayAllocated()  {
    Wrapper* pBlock;
    pBlock = pFirst;
    while (pBlock) {
        printf("%d bytes allocated at %s;%d\n", pBlock->length,
                pBlock->fileName, pBlock->lineNumber);
        pBlock = pBlock->pNext;
    }
}
