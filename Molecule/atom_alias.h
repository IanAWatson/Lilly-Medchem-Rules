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
#ifndef ATOM_ALIAS_H
#define ATOM_ALIAS_H

/*
  From an ISIS reaction file we may have an atom alias 
*/

#include "iwmtypes.h"
#include "iwstring.h"
#include "iwstring_data_source.h"

class Atom_Alias
{
  private:
    atom_number_t _atom;
    IWString      _alias;

//  private functions

    void _copy (const Atom_Alias &);

  public:
    Atom_Alias ();
    Atom_Alias (const Atom_Alias &);

    Atom_Alias & operator= (const Atom_Alias &);

    int build (const const_IWSubstring &, iwstring_data_source &);

    atom_number_t atom_number () const { return _atom;}
    void set_atom_number (atom_number_t s) { _atom = s;}
    const IWString & alias () const { return _alias;}
};

#endif
