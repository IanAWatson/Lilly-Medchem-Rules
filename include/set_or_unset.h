/**************************************************************************

    Copyright (C) 2011  Eli Lilly and Company

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************/
#ifndef IWSET_OR_UNSET_H
#define IWSET_OR_UNSET_H

#include "iwconfig.h"

#include <iostream>

template <typename T>
class Set_or_Unset
{
  protected:
    T _value;
    int _is_set;

  public:
    Set_or_Unset () {_is_set = 0;}
    Set_or_Unset (const T v) {_is_set = 1; _value = v;}

    Set_or_Unset & operator= (const Set_or_Unset &);

    int operator== (const Set_or_Unset &) const;

    int set (const T v) { _value = v; return _is_set = 1;}
    int unset () { return _is_set = 0;}
    int is_set () const { return _is_set;}

    int value (T &) const;
    int matches (const T) const;

    Set_or_Unset<T> & operator = (const T & v) { set (v); return *this;}
};

extern ostream & operator<< (ostream &, const Set_or_Unset<int> &);
extern ostream & operator<< (ostream &, const Set_or_Unset<float> &);
extern ostream & operator<< (ostream &, const Set_or_Unset<double> &);
extern ostream & operator<< (ostream &, const Set_or_Unset<long> &);

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(SET_OR_UNSET_IMPLEMENTATION)

#include <assert.h>

template <typename T>
int
Set_or_Unset<T>::value (T & v) const
{
  if (0 == _is_set)
    return 0;
 
  v = _value;
  return 1;
}

template <typename T>
int
Set_or_Unset<T>::matches (const T v) const
{
  if (! _is_set)
    return 0;

  return _value == v;
}

template <typename T>
ostream &
operator << (ostream & os, const Set_or_Unset<T> & qq)
{
  assert (os.good ());

  T tmp;
  if (qq.value (tmp))
    os << "value is " << tmp;
  else
    os << "value not set.";

  return os;
}

template <typename T>
Set_or_Unset<T> &
Set_or_Unset<T>::operator = (const Set_or_Unset<T> & rhs)
{
  _is_set = rhs._is_set;
  _value = rhs._value;

  return *this;
}

template <typename T>
int
Set_or_Unset<T>::operator == (const Set_or_Unset<T> & rhs) const
{
  if (0 == _is_set && 0 == rhs._is_set)    // both unset
    return 1;

  if (_is_set && rhs._is_set)
    return _value == rhs._value;

  return 0;
}

#endif
#endif
