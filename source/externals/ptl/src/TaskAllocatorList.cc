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
// ---------------------------------------------------------------
//  Tasking class implementation
//

#include <iomanip>

#include "PTL/TaskAllocator.hh"
#include "PTL/TaskAllocatorList.hh"

using namespace PTL;

//======================================================================================//

TaskAllocatorList*&
TaskAllocatorList::fAllocatorList()
{
    static thread_local TaskAllocatorList* _instance = nullptr;
    return _instance;
}

//======================================================================================//

TaskAllocatorList*
TaskAllocatorList::GetAllocatorList()
{
    if(!fAllocatorList())
    {
        fAllocatorList() = new TaskAllocatorList;
    }
    return fAllocatorList();
}

//======================================================================================//

TaskAllocatorList*
TaskAllocatorList::GetAllocatorListIfExist()
{
    return fAllocatorList();
}

//======================================================================================//

TaskAllocatorList::TaskAllocatorList() {}

//======================================================================================//

TaskAllocatorList::~TaskAllocatorList() { fAllocatorList() = nullptr; }

//======================================================================================//

void
TaskAllocatorList::Register(TaskAllocatorBase* alloc)
{
    fList.push_back(alloc);
}

//======================================================================================//

void
TaskAllocatorList::Destroy(int nStat, int verboseLevel)
{
    int    i = 0, j = 0;
    double tmem = 0;
    if(verboseLevel > 0)
    {
        std::cout << "================== Deleting memory pools ==================="
                  << std::endl;
    }
    for(auto& itr : fList)
    {
        double mem = itr->GetAllocatedSize();
        if(i < nStat)
        {
            i++;
            tmem += mem;
            itr->ResetStorage();
            continue;
        }
        j++;
        tmem += mem;
        if(verboseLevel > 1)
        {
            std::cout << "Pool ID '" << itr->GetPoolType()
                      << "', size : " << std::setprecision(3) << mem / 1048576
                      << std::setprecision(6) << " MB" << std::endl;
        }
        itr->ResetStorage();
        // delete *itr;
    }
    if(verboseLevel > 0)
    {
        std::cout << "Number of memory pools allocated: " << Size()
                  << "; of which, static: " << i << std::endl;
        std::cout << "Dynamic pools deleted: " << j
                  << " / Total memory freed: " << std::setprecision(2) << tmem / 1048576
                  << std::setprecision(6) << " MB" << std::endl;
        std::cout << "============================================================"
                  << std::endl;
    }
    fList.clear();
}

//======================================================================================//

int
TaskAllocatorList::Size() const
{
    return fList.size();
}

//======================================================================================//
