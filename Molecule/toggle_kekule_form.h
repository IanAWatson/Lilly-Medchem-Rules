/**************************************************************************

    Copyright (C) 2012  Eli Lilly and Company

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
#ifndef TOGGLE_KEKULE_FORM_H
#define TOGGLE_KEKULE_FORM_H

#include "molecule.h"

class Command_Line;
class msi_attribute;
class Bond;

/*
  We want to fix a number of molecules into a single Kekule form
*/

class Toggle_Kekule_Form_Temporary_Arrays
{
  private:
//  At the end of the computation we recompute the aromaticity to see if
//  all our atoms are still aromatic

    aromaticity_type_t * _arom;

//  For efficiency, we get a copy of the atoms in the molecule

    Atom ** _atom;

//  We keep track of whether or not each atom already has a double bond

    int * _has_double_bond;

//  We also keep track of which bonds can change
  
    int * _can_change_bond;

    int * _process_these;

//  When we get an atom that cannot be changed, all rings that contain
//  that atom cannot change

    int * _ring_can_vary;

    int * _atom_can_change;

  public:
    Toggle_Kekule_Form_Temporary_Arrays (Molecule &);
    ~Toggle_Kekule_Form_Temporary_Arrays ();

    int * can_change_bond () { return _can_change_bond;}
    aromaticity_type_t * arom () { return _arom;}
    Atom ** atom () { return _atom;}
    int * has_double_bond () { return _has_double_bond;}
    int * process_these () { return _process_these;}

    void set_ring_can_toggle (int r, int s) { _ring_can_vary[r] = s;}
    int  ring_can_toggle (int r) const { return _ring_can_vary[r];}

    void set_atom_can_change (int a, int s) { _atom_can_change[a] = s;}
    int  atom_can_change (int a) const { return _atom_can_change[a];}
};

class Toggle_Kekule_Form
{
  private:
    resizable_array_p<Bond> _bond;

//  Whether or not each Bond is already of the requested bond type

    int * _correct;

//  By default, we do NOT allow a pyrrole Nitrogen to toggle

    int _allow_pyrrole_to_change;

    int _display_error_messages;

//  private functions

    void _no_changes_to_atom (Molecule & m,
                             atom_number_t zatom,
                             Toggle_Kekule_Form_Temporary_Arrays &) const;

    void _set_all_bonds_to_single (Molecule & m,
                                   int id,
                                   Toggle_Kekule_Form_Temporary_Arrays &) const;

    int _bond_is_correct (const Molecule & m,
                          const Set_of_Atoms & embedding,
                          const Bond * b) const;
    int _all_bonds_correct (const Molecule & m,
                            const Set_of_Atoms & embedding);
    int _all_bonds_aromatic (Molecule & m,
                             const Set_of_Atoms & embedding) const;
    int _all_atoms_aromatic (Molecule & m, int, int, Toggle_Kekule_Form_Temporary_Arrays &) const;
    int _ring_is_involved (const Ring * r) const;

    int _get_ring_system_atoms (resizable_array<int> & atoms_to_process,
                       int rid,
                       atom_number_t zatom,
                       Toggle_Kekule_Form_Temporary_Arrays & ) const;
    int _set_our_bonds (Molecule & m,
                        const Set_of_Atoms & embedding,
                        int id,
                        Toggle_Kekule_Form_Temporary_Arrays &) const;

    int _do_not_process_rings_containing (Molecule & m,
                                         atom_number_t zatom,
                                         Toggle_Kekule_Form_Temporary_Arrays & tkfta) const;

//  void _do_chemistry (Molecule & m,
//                      int id,
//                      Toggle_Kekule_Form_Temporary_Arrays &) const;
    void _do_chemistry (Molecule & m,
                        Toggle_Kekule_Form_Temporary_Arrays &) const;
    void _do_chemistry (Molecule & m,
                        const Ring & r,
                        Toggle_Kekule_Form_Temporary_Arrays &) const;
    void _do_chemistry_aromatic_ring (Molecule & m, const Ring & r, Toggle_Kekule_Form_Temporary_Arrays & tkfta) const;
    int _process_ring_system (Molecule & m,
                              resizable_array<int> & atoms_to_process,
                              int rid,
                              int zitem,
                              Toggle_Kekule_Form_Temporary_Arrays &) const;
    int _process_ring_system (Molecule & m,
                              const Set_of_Atoms & embedding,
                              int rid,
                              Toggle_Kekule_Form_Temporary_Arrays &);
    int _process (Molecule & m,
                  const int * process_these,
                  int * already_done);
    int _process (Molecule & m, const Set_of_Atoms & embedding, int * process_these);
    int _process (Molecule & m,
                  const Set_of_Atoms & embedding,
                  Toggle_Kekule_Form_Temporary_Arrays &);
    int _process_single_ring (Molecule & m,
                              const Set_of_Atoms & embedding,
                              int id,
                              Toggle_Kekule_Form_Temporary_Arrays &);
    int _process_single_ring2 (Molecule & m,
                              const Set_of_Atoms & embedding,
                              int id,
                              Toggle_Kekule_Form_Temporary_Arrays &);

  public:
    Toggle_Kekule_Form ();
    ~Toggle_Kekule_Form ();

    int ok () const;
    int debug_print (ostream &) const;

    int active () const { return _bond.number_elements ();}

    void set_display_error_messages (int s) { _display_error_messages = s;}

    const Bond * contains_bond (atom_number_t a1, atom_number_t a2) const;
    int will_change_ring (const Ring * r, const Set_of_Atoms &) const;

    int construct_from_command_line (Command_Line &, char, int = 0);
    int add_bond_from_msi_attribute (const msi_attribute &);
    int add_bond (int, int, bond_type_t);
    int add_bond (Bond * b);

    void set_allow_pyrrole_to_change (int s) { _allow_pyrrole_to_change = s;}

    int write_msi (ostream &, const IWString &, const char *) const;

    int ok_embedding (const Set_of_Atoms & embedding) const;

    int process (Molecule &, const Set_of_Atoms &, int &);

    int process (Molecule &, atom_number_t a1, atom_number_t a2, bond_type_t bt, int & changed);
};

#endif
