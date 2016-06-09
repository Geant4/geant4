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
void tree<N>::iterator::clear(void)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  for (sibling_iterator i(beginChildren());
       i != endChildren(); )
  {
    i.clear();
    ++i;
    (i - 1).destroy();
  }

  current->first_child = NULL;
  current->last_child = NULL;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::erase(void)
{
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  /* remove all children and grandchildren */
  clear();

  /* make copies because there is an ordering
     problem updating the previous and next links */
  tree_node<N>* old_prev = current->prev_sibling;
  tree_node<N>* old_next = current->next_sibling;

  if (old_prev != NULL)
    old_prev->next_sibling = old_next;

  if (old_next != NULL)
    old_next->prev_sibling = old_prev;

  if (current->parent->first_child == current)
    current->parent->first_child = old_next;

  if (current->parent->last_child == current)
    current->parent->last_child = old_prev;

  /* destroy has no parameters, so gcc complains that it
     can't find the method unless this is made explicit */
  super::destroy();

  return iterator(old_next);
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::insert(const N& data)
{
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  tree_node<N>* p = new tree_node<N>(data);
  p->parent = current->parent;
  p->first_child = p->last_child = NULL;

  p->prev_sibling = current->prev_sibling;
  p->next_sibling = current;

  if (current->prev_sibling != NULL)
    current->prev_sibling->next_sibling = p;
  else
    current->parent->first_child = p;

  current->prev_sibling = p;

  return p;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::insert(const tree<N>& t)
{
  assert(this->current != NULL);

  return absorb_before(new tree<N>(t));
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::insert_after(const N& data)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  tree_node<N>* p = new tree_node<N>(data);
  p->parent = current->parent;
  p->first_child = p->last_child = NULL;

  p->prev_sibling = current;
  p->next_sibling = current->next_sibling;

  if (current->next_sibling != NULL)
    current->next_sibling->prev_sibling = p;
  else
    current->parent->last_child = p;

  current->next_sibling = p;

  return p;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::insert_after(const tree<N>& t)
{
  assert(this->current != NULL);

  return absorb_after(new tree<N>(t));
}

template<class N>
void tree<N>::iterator::pop_back(void)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);
  assert(current->last_child != NULL);
  assert(current->last_child->next_sibling == NULL);
  
  tree_node<N>* p = current->last_child;
  current->last_child = p->prev_sibling;

  if (p->prev_sibling != NULL)
    p->prev_sibling->next_sibling = NULL;
  else
    current->first_child = NULL;

  delete p;
}

template<class N>
void tree<N>::iterator::pop_front(void)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);
  assert(current->first_child != NULL);
  assert(current->first_child->prev_sibling == NULL);

  tree_node<N>* p = current->first_child;
  current->first_child = p->next_sibling;

  if (p->next_sibling != NULL)
    p->next_sibling->prev_sibling = NULL;
  else
    current->last_child = NULL;

  delete p;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::push_back(const N& data)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  tree_node<N>* p = new tree_node<N>(data);
  p->parent = current;
  p->first_child = p->last_child = NULL;

  if (current->last_child == NULL)
  {
    current->first_child = current->last_child = p;
    p->prev_sibling = p->next_sibling = NULL;
  }
  else
  {
    p->prev_sibling = current->last_child;
    p->next_sibling = NULL;
    current->last_child->next_sibling = p;
    current->last_child = p;
  }

  return iterator(p);
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::push_back(const tree<N>& t)
{
  assert(this->current != NULL);

  return absorb_back(new tree<N>(t));
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::push_front(const N& data)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  tree_node<N>* p = new tree_node<N>(data);
  p->parent = current;
  p->first_child = p->last_child = NULL;

  if (current->first_child == NULL)
  {
    current->first_child = current->last_child = p;
    p->prev_sibling = p->next_sibling = NULL;
  }
  else
  {
    p->prev_sibling = NULL;
    p->next_sibling = current->first_child;
    current->first_child->prev_sibling = p;
    current->first_child = p; 
  }

  return iterator(p);
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::push_front(const tree<N>& t)
{
  assert(this->current != NULL);

  return absorb_front(new tree<N>(t));
}

template<class N>
N tree<N>::iterator::replace(const N& data)
{
  N ret(this->current->data);
  this->current->data = data;
  return ret;
}

template<class N>
size_t tree<N>::iterator::size(void) const
{
  assert(this->current != NULL);

  size_t s = 1;

  for (typename tree<N>::sibling_iterator i(beginChildren());
       i != endChildren(); ++i)
  {
    s += i.size();
  }

  return s;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::splice_after(tree<N>& t)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);

  tree_node<N>* p = t.master->next_sibling;
  t.master->next_sibling = t.master;
  t.master->prev_sibling = t.master;

  p->parent = current->parent;

  p->prev_sibling = current;
  p->next_sibling = current->next_sibling;

  if (current->next_sibling != NULL)
    current->next_sibling->prev_sibling = p;
  else
    current->parent->last_child = p;

  current->next_sibling = p;

  return p;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::splice_back(tree<N>& t)
{
  /* this line necessary for g++ 3.4 onward */
  tree_node<N>*& current = this->current;

  assert(current != NULL);
  assert(!t.empty());

  tree_node<N>* p = t.master->next_sibling;
  t.master->next_sibling = t.master;
  t.master->prev_sibling = t.master;

  p->parent = current;

  if (current->last_child == NULL)
  {
    current->first_child = current->last_child = p;
    p->prev_sibling = p->next_sibling = NULL;
  }
  else
  {
    p->prev_sibling = current->last_child;
    p->next_sibling = NULL;
    current->last_child->next_sibling = p;
    current->last_child = p;
  }

  return p;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::splice_before(tree<N>& t)
{
  /* TODO */
  assert(false);
  return NULL;
}

template<class N>
typename tree<N>::iterator tree<N>::iterator::splice_front(tree<N>& t)
{
  /* TODO */
  assert(false);
  return NULL;
}

