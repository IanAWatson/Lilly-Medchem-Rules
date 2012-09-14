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
#ifndef MOLECULE_SMARTS_H
#define MOLECULE_SMARTS_H

  private:
    int _smarts (atom_number_t astart,
                   int * include_atom,
                   int flag,
                   IWString & s);

    void _compute_ncon_and_explicit_hydrogens (atom_number_t zatom,
                                                int & ncon,
                                                int & eh,
                                                const int * include_atom) const;
    void _append_isotope_and_atomic_symbol (atom_number_t zatom,
                                             IWString & smiles);
    int _append_smarts_equivalent_for_atom (atom_number_t zatom,
                                              int ncon,
                                              int rm,
                                              IWString & s) const;


#endif

/* arch-tag: 12f3775a-c3d9-4d19-a884-b63de1606b98

*/
