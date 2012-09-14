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
    int _parse_tripos_atom_record (const const_IWSubstring & buffer, atom_number_t, int &);
    int _parse_tripos_bond_record (const const_IWSubstring & buffer, int * aromatic_atoms, int, int * aromatic_bond);
    int _mol2_assign_default_formal_charges ();
    int _doubly_bonded_to_oxygen (atom_number_t zatom) const;
    int _tripos_atom_type_from_string (atom_number_t, const const_IWSubstring &);
    int _place_formal_charges_on_quat_n_from_mol2 ();
