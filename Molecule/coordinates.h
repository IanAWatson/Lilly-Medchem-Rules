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
#ifndef ATOM_COORDINATES_H
#define ATOM_COORDINATES_H

#include "space_vector.h"
#include "iwmtypes.h"

class Atom;

class Coordinates : public Space_Vector<coord_t>
{
  private:

  public:
    Coordinates () {};
    Coordinates (const Space_Vector<coord_t> &);

    Coordinates (coord_t cx, coord_t cy, coord_t cz) : Space_Vector<coord_t> (cx, cy, cz) {};
    Coordinates (const Coordinates & r) : Space_Vector<coord_t> (r.x (), r.y (), r.z ()) {};
    Coordinates (const Atom & a);
};

class Coordinates_double : public Space_Vector<double>
{
  private:

  public:
    Coordinates_double () {};
    Coordinates_double (const Space_Vector<double> &);

    Coordinates_double (double cx, double cy, double cz) : Space_Vector<double> (cx, cy, cz) {};
    Coordinates_double (const Coordinates & r) : Space_Vector<double> (r.x (), r.y (), r.z ()) {};
    Coordinates_double (const Atom & a);
};

#endif
