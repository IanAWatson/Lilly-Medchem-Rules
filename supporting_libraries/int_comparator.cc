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

#include "misc.h"

int
Int_Comparator_Larger::operator() (int i1, int i2) const
{
  if (i1 < i2)
    return -1;

  if (i1 > i2)
    return 1;

  return 0;
}

int
Int_Comparator_Smaller::operator() (int i1, int i2) const
{
  if (i1 < i2)
    return 1;

  if (i1 > i2)
    return -1;

  return 0;
}
