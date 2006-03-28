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

template<class N>
template<class T, class R, class P>
void tree<N>::tree_iterator<T, R, P>::post_inc(tree_node<N>*& current, long& depth)
{
  assert(current != NULL);

  if (current->next_sibling == NULL)
  {
    current = current->parent;
    depth--;
  }
  else
  {
    current = current->next_sibling;
    while (current->first_child != NULL)
    {
      current = current->first_child;
      depth++;
    }
  }
}

template<class N>
template<class T, class R, class P>
void tree<N>::tree_iterator<T, R, P>::post_dec(tree_node<N>*& current, long& depth)
{
  assert(current != NULL);

  if(current->last_child == NULL)
  {
    while (current->prev_sibling == NULL)
    {
      current = current->parent;
      depth--;
    }

    current = current->prev_sibling;
  }
  else
  {
    current = current->last_child;
    depth++;
  }
}

template<class N>
template<class T, class R, class P>
void tree<N>::tree_iterator<T, R, P>::pre_inc(tree_node<N>*& current, long& depth)
{
  assert(current != NULL);

  if (current->first_child != NULL)
  {
    current = current->first_child;
    depth++;
  }
  else
  {
    while (current->next_sibling == NULL)
    {
      current = current->parent;
      depth--;
      if (current == NULL)
        return;
    }
    current = current->next_sibling; 
  }
}

template<class N>
template<class T, class R, class P>
void tree<N>::tree_iterator<T, R, P>::pre_dec(tree_node<N>*& current, long& depth)
{
  assert(current != NULL);

  if (current->prev_sibling)
  {
    current = current->prev_sibling;
    while (current->last_child != NULL)
    {
      current = current->last_child;
      depth++;
    }
  }
  else
  {
    current = current->parent;
    depth--;
  }
}

template<class N>
template<class T, class R, class P>
inline void tree<N>::tree_iterator<T, R, P>::sib_inc(tree_node<N>*& current)
{
  if (current != NULL)
    current = current->next_sibling;
}

template<class N>
template<class T, class R, class P>
inline void tree<N>::tree_iterator<T, R, P>::sib_dec(tree_node<N>*& current)
{
  if (current != NULL)
    current = current->prev_sibling;
}
