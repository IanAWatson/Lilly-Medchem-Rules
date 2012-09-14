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
/*
  Some extra private functions needed by careful_frag
*/

#ifndef CAREFUL_FRAG_H
#define CAREFUL_FRAG_H

    int _reduce_to_largest_fragment_carefully (Fragment_Data * fc, int * already_counted);
    int _is_nitro (atom_number_t, int *) const;
    int _is_sulphate_like (atom_number_t, int *) const;
    int _identify_fragment_undesirable_groups (int * exclude) const;
#endif
