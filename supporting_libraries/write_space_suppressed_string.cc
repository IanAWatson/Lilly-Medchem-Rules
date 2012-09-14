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
#include "misc.h"
#include "iwstring.h"

ostream &
write_space_suppressed_string (const IWString & zstring,
                               ostream & os,
                               char fill_char)
{
  if (zstring.contains (' '))
  {
    IWString tmp (zstring);
    tmp.gsub (' ', fill_char);

    os << tmp;
  }
  else
    os << zstring;

  return os;
}

