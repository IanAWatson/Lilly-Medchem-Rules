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
#ifndef MDL_ATOM_RECORD_H
#define MDL_ATOM_RECORD_H

#include "molecule.h"

class MDL_Atom_Record
{
  private:
    coord_t _x;
    coord_t _y;
    coord_t _z;
    const_IWSubstring _atomic_symbol;
    int _msdiff;
    int _chg;
    int _astere;
    int _hcount;
    int _stereo_care;
    int _valence;
    int _h0designator;
    int _atom_map;
    int _inversion;
    int _exact_change;

//  private functions

    void _default_values();

  public:
    MDL_Atom_Record();

    int build (const const_IWSubstring &);

    Atom * create_atom() const;

    int astere() const { return _astere;}
    int hcount() const { return _hcount;}
    int valence() const { return _valence;}
    int h0designator() const { return _h0designator;}

    int atom_map() const { return _atom_map;}
    int inversion() const { return _inversion;}
    int exact_change() const { return _exact_change;}

    const const_IWSubstring & atomic_symbol() const { return _atomic_symbol;}
};

class MDL_Bond_Record
{
  private:
    atom_number_t _a1;
    atom_number_t _a2;
    int _bond_type_read_in;
    int _directionality;
    int _bond_stereo;
    int _bond_topology;
    int _reacting_center_status;

  public:
    MDL_Bond_Record();

    int build (const const_IWSubstring &, int natoms);

    atom_number_t a1() const { return _a1;}
    atom_number_t a2() const { return _a2;}

    int bond_type_read_in() const { return _bond_type_read_in;}

    int directionality() const { return _directionality;}

    int bond_topology() const { return _bond_topology;}

    int reacting_center_status() const { return _reacting_center_status;}

    int convert_mdl_bond_type_read_in_to_query (bond_type_t & for_query,
                                        bond_type_t & for_building_a_molecule) const;

//  Based on the bond type in the input file, we produce one bond to be
//  used for building a molecule, and another bond type that encodes
//  all the AND and OR conditions

    int bond_type_for_molecule(bond_type_t &) const;
    int bond_type_for_query(bond_type_t &) const;
};

#endif
