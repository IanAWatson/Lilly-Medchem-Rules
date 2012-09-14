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
#ifndef COLLECTION_TEMPLATE_H
#define COLLECTION_TEMPLATE_H

/*
  There are a couple of instances where we need a collection of
  objects as well as a text description of what the items are.
  The first two examples are partial charges, and atom types
*/

#include "iwstring.h"

template <typename T>
class Collection_Template : public resizable_array<T>
{
  private:
    IWString _type;

  public:

    Collection_Template<T> & operator = (const Collection_Template<T> & rhs)
      {
        resizable_array<T>::operator= (rhs);
        _type = rhs._type;

        return *this;
      }

    IWString & ztype () { return _type;}
    const IWString & ztype () const { return _type;}

    void set_type (const IWString & t) { _type = t;}
};

#endif
