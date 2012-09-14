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

#include "accumulator.h"


KahanSum::KahanSum ()
{
  _c = 0.0;
  _sum = 0.0;

  return;
}

KahanSum &
KahanSum::operator = (double d)
{
  _sum = d;
  _c = 0.0;

  return *this;
}

KahanSum &
KahanSum::operator += (double d)
{
  double y = d - _c;
  double t = _sum + y;
  _c = (t - _sum) - y;
  _sum = t;

  return *this;
}

KahanSum &
KahanSum::operator += (const KahanSum & rhs)
{
  _c += rhs._c;

  return this->operator+=(rhs._sum);
}
