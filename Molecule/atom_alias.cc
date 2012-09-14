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
#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

#include "atom_alias.h"

Atom_Alias::Atom_Alias ()
{
  _atom = INVALID_ATOM_NUMBER;

  return;
}

Atom_Alias::Atom_Alias (const Atom_Alias & rhs)
{
  _copy (rhs);

  return;
}

Atom_Alias &
Atom_Alias::operator= (const Atom_Alias & rhs)
{
  _copy (rhs);

  return *this;
}

void
Atom_Alias::_copy (const Atom_Alias & rhs)
{
  _atom = rhs._atom;
  _alias = rhs._alias;

  return;
}

int
Atom_Alias::build (const const_IWSubstring & buffer,
                   iwstring_data_source & input)
{
  assert (buffer.starts_with ("A  "));

  const_IWSubstring tmp (buffer);

  tmp.remove_leading_chars (3);

  tmp.strip_leading_blanks ();

  if (! tmp.numeric_value (_atom) || _atom < 1)
  {
    cerr << "Atom_Alias::build: invalid atom number specification '" << buffer << "'\n";
    return 0;
  }

  _atom--;

  if (! input.next_record (_alias))
  {
    cerr << "Atom_Alias::build:premature EOF\n";
    return 0;
  }

  return 1;
}

template class resizable_array_p<Atom_Alias>;
template class resizable_array_base<Atom_Alias *>;

