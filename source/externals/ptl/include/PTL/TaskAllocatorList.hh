//
// MIT License
// Copyright (c) 2020 Jonathan R. Madsen
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
// "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// ------------------------------------------------------------
// Tasking class header file
//
// Class Description:
//
// A class to store all TaskAllocator objects in a thread for the sake
// of cleanly deleting them.
//
// ------------------------------------------------------------

#pragma once

#include "PTL/Globals.hh"
#include <vector>

namespace PTL
{
class TaskAllocatorBase;

class TaskAllocatorList
{
public:  // with description
    static TaskAllocatorList* GetAllocatorList();
    static TaskAllocatorList* GetAllocatorListIfExist();

public:
    ~TaskAllocatorList();
    void Register(TaskAllocatorBase*);
    void Destroy(int nStat = 0, int verboseLevel = 0);
    int  Size() const;

private:
    TaskAllocatorList();

private:
    static TaskAllocatorList*&      fAllocatorList();
    std::vector<TaskAllocatorBase*> fList;
};

}  // namespace PTL
