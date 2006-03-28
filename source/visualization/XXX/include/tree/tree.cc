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

template <class N>
template <class InputIterator>
tree<N>::tree(InputIterator first, InputIterator last, unsigned int n)
{
  assert(n != 0);

  createMaster();

  if (first == last)
    return;

  setRoot(*first);

  /* don't directly initialize i to first + 1 because not all STL
     iterators support addition */
  InputIterator i(first);
  
  if (++i == last)
    return;

  /* Make a list of leaves. Initially the list contains only the
     root.  Pick a leaf and add items to it until it is full.
     Added items become new leaves.  Continue until there
     are no more items */

  std::list<typename tree<N>::iterator> leaves;
  leaves.push_back(getRoot());

  for (;;)
  {
    typename tree<N>::iterator leaf(leaves.front());

    for (unsigned int k = 0; k < n; ++k)
    {
      leaves.push_back(leaf.push_back(*i));

      if (++i == last)
        return; 
    }

    leaves.pop_front();
  }
}

template <class N>
void tree<N>::copy_subtree(tree<N>::const_iterator subtree)
{
  /* Copying replaces this tree. */
  clear();

  /* Original is empty, so leave it empty. */
  if (subtree == NULL)
    return;

  /* Copy the root data. */
  sibling_iterator root(setRoot(*subtree));

  /* Recursively copy the children. */ 
  for (const_sibling_iterator i(subtree.beginChildren());
       i != subtree.endChildren(); ++i)
    root.absorb_back(new tree<N>(i));
}

template <class N>
void tree<N>::clear(void)
{
  if (empty())
    return;

  /* delete all nodes bottom-up */
  for (postorder_iterator i(beginPost()); i != endPost(); )
  {
    postorder_iterator tmp(i);
    ++i;
    tmp.destroy();
  }

  /* Preserve master invariant. */
  master->first_child  = master->last_child   = NULL;
  master->prev_sibling = master->next_sibling = master; 

  /* Leave the master for the destructor. */
}

template <class N>
size_t tree<N>::height(void) const
{
  size_t h = 0;

  for (const_preorder_iterator i(beginPre()); i != endPre(); ++i)
  {
    if (static_cast<size_t>(i.depth() + 1) > h)
      h = i.depth() + 1;
  }

  return h;
}

template <class N>
size_t tree<N>::height_leftmost(void) const
{
  size_t h = 0;

  tree_node<N>* p = master->next_sibling;

  if (p != master)
  {
    do {
      ++h;
      p = p->first_child;
    } while (p != NULL);
  }  

  return h;
}

template<class N>
size_t tree<N>::size(void) const
{
  size_t s = 0;

  for (const_preorder_iterator i(beginPre()); i != endPre(); ++i)
    ++s;  

  return s;
} 
