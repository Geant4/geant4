/*
Copyright (c) 2003-2006, Troy Aaron Johnson
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of Troy Aaron Johnson nor the names of any
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <string>

#include "tree.h"

using namespace taj;

bool pruning_predicate(const tree<std::string>::iterator& node)
{
  return (*node == "B");
}

template <class TreeType>
void print(const tree<TreeType>& t)
{
  std::cout << "preorder traversal: " << std::endl;
  for (typename tree<TreeType>::const_preorder_iterator i = t.beginPre(); i != t.endPre(); ++i)
  {
    std::cout << *i << std::endl;
  }
  std::cout << std::endl;

  std::cout << "postorder traversal: " << std::endl;
  for (typename tree<TreeType>::const_postorder_iterator i = t.beginPost(); i != t.endPost(); ++i)
  {
    std::cout << *i << std::endl;
  }
  std::cout << std::endl;

  std::cout << "bfs traversal: " << std::endl;
  for (typename tree<TreeType>::const_bfs_iterator i = t.beginBfs(); i != i.endBfs(); ++i)
  {
    std::cout << *i << std::endl;
  }
  std::cout << std::endl;
}

void pointer_specialization_test(void)
{
  tree<int*> pint_tree;

  int x = 1, y = 2, z = 3;

  tree<int*>::iterator i(pint_tree.setRoot(&x));

  i.push_back(&y);
  i.push_back(&z);

  std::cout << "int* tree root is " << **i << std::endl;

  tree<int*>::sibling_iterator si(i.beginChildren());

  tree<char*> t1;
  tree<long*> t2;
  tree<float*> t3;
  tree<double*> t4;
  tree<std::string*> t5;
  tree<unsigned char*> t6;
  tree<unsigned long*> t7;
}

int main(int argc, char* argv[])
{
  tree<std::string> t;

  assert(t.empty());
  assert(t.size() == 0);

  tree<std::string>::sibling_iterator i(t.setRoot("root"));

  assert(!t.empty());
  assert(t.size() == 1);
  assert(*i == "root");

  *i = "new root";

  assert(*i == "new root");

  tree<std::string>::sibling_iterator j(i.push_back("left"));
  assert(t.size() == 2);
  tree<std::string>::sibling_iterator k(i.push_back("right"));
  assert(t.size() == 3);

  assert(*j == "left");
  assert(*k == "right");

  j.push_back("A");
  j.push_back("B");

  k.push_back("C");
  k.push_back("D");

  assert(find(t.beginPre(), t.endPre(), "C") != t.endPre());
  assert(find(t.beginPre(), t.endPre(), "Z") == t.endPre());

  std::cout << "tree height is " << t.height() << std::endl;

  print(t);

  tree<std::string> copy_of_t(t);

  print(copy_of_t);

  std::list<int> int_list;
  for (int i = 1; i < 20; ++i)
    int_list.push_back(i);

  tree<int> int_tree(int_list.begin(), int_list.end(), 2);

  print(int_tree);

  pointer_specialization_test();
}
