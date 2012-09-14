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
#include <stdlib.h>
#include <memory>
#include <limits>
#include <ctype.h>

using namespace std;

#include "cmdline.h"
#include "misc.h"
#include "iw_auto_array.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "msi_object.h"

#include "toggle_kekule_form.h"
#include "path.h"

Toggle_Kekule_Form::Toggle_Kekule_Form ()
{
  _correct = NULL;

  _allow_pyrrole_to_change = 0;

  _display_error_messages = 1;

  return;
}

Toggle_Kekule_Form::~Toggle_Kekule_Form ()
{
  if (NULL != _correct)
    delete [] _correct;

  return;
}

static int
fetch_int (const const_IWSubstring & c, int & zresult)
{
  zresult = 0;

  for (int i = 0; i < c.nchars (); i++)
  {
    if (! isdigit (c[i]))
      return i;

    zresult = 10 * zresult + c[i] - '0';
  }

  return c.nchars ();
}

int
Toggle_Kekule_Form::debug_print (ostream & os) const
{
  cerr << "Toggle_Kekule_Form with " << _bond.number_elements () << " bonds\n";
  for (int i = 0; i < _bond.number_elements (); i++)
  {
    const Bond * b = _bond[i];

    cerr << "Between matched atoms " << b->a1 () << " and " << b->a2 () << " type ";
    if (b->is_single_bond ())
      cerr << "single";
    else if (b->is_double_bond ())
      cerr << "double";
    else
      cerr << "what kind of bond is this " << *b;

    cerr << endl;
  }

  return os.good ();
}

const Bond *
Toggle_Kekule_Form::contains_bond (atom_number_t a1, atom_number_t a2) const
{
  for (int i = 0; i < _bond.number_elements (); i++)
  {
    if (_bond[i]->involves (a1, a2))
      return _bond[i];
  }

  return NULL;
}

/*
  Do we have any bonds which are in a given ring
*/

int
Toggle_Kekule_Form::will_change_ring (const Ring * r,
                                      const Set_of_Atoms & embedding) const
{
  int nb = _bond.number_elements ();

  for (int i = 0; i < nb; i++)
  {
    int i1 = _bond[i]->a1 ();
    int i2 = _bond[i]->a2 ();

    atom_number_t a1 = embedding[i1];
    atom_number_t a2 = embedding[i2];

    if (r->contains_bond (a1, a2))
      return 1;
  }
  
  return 0;     // none of our bonds are in this ring
}

int
Toggle_Kekule_Form::construct_from_command_line (Command_Line & cl, 
                                                 char c,
                                                 int verbose)
{
  const_IWSubstring cv;
  int i = 0;
  while (cl.value (c, cv, i++))
  {
    int a1;
    int nchars = fetch_int (cv, a1);
    if (0 == nchars || cv.nchars () == nchars)
    {
      cerr << "Toggle_Kekule_Form::construct_from_command_line: invalid a1 '" << cv << "'\n";
      return 0;
    }

    bond_type_t bt;
    if ('-' == cv[nchars])
      bt = SINGLE_BOND;
    else if ('=' == cv[nchars])
      bt = DOUBLE_BOND;
    else
    {
      cerr << "Toggle_Kekule_Form::construct_from_command_line: bonds must be '-' or '=', '" << cv << "' is invalid\n";
      return 0;
    }

    cv += nchars + 1;

    int a2;
    nchars = fetch_int (cv, a2);
    if (0 == nchars || nchars != cv.nchars ())
    {
      cerr << "Toggle_Kekule_Form::construct_from_command_line: invalid a1 '" << cv << "'\n";
      return 0;
    }

    if (a1 == a2)
    {
      cerr << "Toggle_Kekule_Form::construct_from_command_line:atoms must be distinct\n";
      return 0;
    }

    Bond * b = new Bond (a1, a2, bt);

    _bond.add (b);
  }

  _correct = new int[_bond.number_elements ()];

  return 1;
}

int
Toggle_Kekule_Form::ok_embedding (const Set_of_Atoms & embedding) const
{
  int rc = 1;
  for (int i = 0; i < _bond.number_elements (); i++)
  {
    const Bond * b = _bond[i];
    if (b->a1 () >= embedding.number_elements ())
    {
      cerr << "Toggle_Kekule_Form::ok_embedding: bond " << i << " atom " << b->a1 () << " invalid\n";
      rc = 0;
    }

    if (b->a2 () >= embedding.number_elements ())
    {
      cerr << "Toggle_Kekule_Form::ok_embedding: bond " << i << " atom " << b->a2 () << " invalid\n";
      rc = 0;
    }
  }

  if (0 == rc)
    cerr << embedding << endl;

  return rc;
}

int
Toggle_Kekule_Form::_bond_is_correct (const Molecule & m,
                                      const Set_of_Atoms & embedding,
                                      const Bond * b) const
{
  atom_number_t a1 = embedding[b->a1 ()];
  atom_number_t a2 = embedding[b->a2 ()];

  const Bond * existing_bond = m.bond_between_atoms (a1, a2);
  assert (NULL != existing_bond);

//#define DEBUG_BOND_IS_CORRECT
#ifdef DEBUG_BOND_IS_CORRECT
  cerr << "Matched atoms " << b->a1 () << " and " << b->a2 () << ", atoms " << a1 << " and " << a2 << endl;
  cerr << "Checking bond    " << (*b) << endl;
  cerr << "against existing " << (*existing_bond) << endl;
#endif

  if (b->is_single_bond () && existing_bond->is_single_bond ())
    return 1;

  if (b->is_double_bond () && existing_bond->is_double_bond ())
    return 1;

  return 0;
}

int
Toggle_Kekule_Form::_all_bonds_correct (const Molecule & m,
                                        const Set_of_Atoms & embedding)
{
  int nb = _bond.number_elements ();

  int rc = 1;

  for (int i = 0; i < nb; i++)
  {
    if (_bond_is_correct (m, embedding, _bond[i]))
      _correct[i] = 1;
    else
    {
      _correct[i] = 0;
      rc = 0;
    }
  }

#ifdef DEBUG_BOND_IS_CORRECT
  cerr << "Toggle_Kekule_Form::_all_bonds_correct:returning " << rc << endl;
#endif

  return rc;
}

int
Toggle_Kekule_Form::_all_bonds_aromatic (Molecule & m,
                                        const Set_of_Atoms & embedding) const
{
  int nb = _bond.number_elements ();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond[i];
    atom_number_t a1 = embedding[b->a1 ()];
    atom_number_t a2 = embedding[b->a2 ()];

    if (! m.in_same_aromatic_ring (a1, a2))
    {
      cerr << "Matched atoms " << a1 << " and " << a2 << " not in the same aromatic ring\n";
      cerr << "Bond " << i << " item " << b->a1 () << " and item " << b->a2 () << " in " << embedding << endl;
      return 0;
    }
  }

  return 1;
}

int
Toggle_Kekule_Form::_ring_is_involved (const Ring * r) const
{
  int nb = _bond.number_elements ();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond[i];

    if (r->contains_bond (b->a1 (), b->a2 ()))
      return 1;
  }

  return 0;
}

int
Toggle_Kekule_Form::_do_not_process_rings_containing (Molecule & m,
                                         atom_number_t zatom,
                                         Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  tkfta.set_atom_can_change(zatom, 0);

  int * can_change_bond = tkfta.can_change_bond ();

  const int nr = m.nrings ();

  const int matoms = m.natoms();

  for (int i = 0; i < nr; i++)
  {
    if (! tkfta.ring_can_toggle(i))
      continue;

    const Ring * ri = m.ringi (i);

    if (! ri->contains (zatom))
      continue;

    int ring_size = ri->number_elements();

    tkfta.set_ring_can_toggle (i, 0);

    atom_number_t prev = ri->item(0);

    tkfta.set_atom_can_change(prev, 0);

    for (int j = 1; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      tkfta.set_atom_can_change(k, 0);

      can_change_bond[prev * matoms + k] = can_change_bond[k * matoms + prev] = 0;

      if (2 == m.nrings(prev) && 2 == m.nrings(k))
      {
        _do_not_process_rings_containing(m, prev, tkfta);
        _do_not_process_rings_containing(m, k, tkfta);
      }

      prev = k;
    }
  }

  return 1;
}

/*
  Mark all the bonds to atom ZATOM as non-changing
*/

void
Toggle_Kekule_Form::_no_changes_to_atom (Molecule & m,
                                         atom_number_t zatom,
                                         Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
//#define DEBUG_NO_CHANGES_TO_ATOM
#ifdef DEBUG_NO_CHANGES_TO_ATOM
  cerr << "Atom " << zatom << " '" << m.smarts_equivalent_for_atom (zatom) << "' cannot change\n";
#endif

  tkfta.set_atom_can_change (zatom, 0);

  int * can_change_bond = tkfta.can_change_bond ();
  Atom * const * atom = tkfta.atom ();

  int matoms = m.natoms ();

  const Atom * a = atom[zatom];

  int acon = a->ncon ();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (zatom, i);

    can_change_bond[zatom * matoms + j] = can_change_bond[j * matoms + zatom] = 0;
//  cerr << "No changes to bond between " << zatom << " and " << j << endl;
  }

  tkfta.set_atom_can_change (zatom, 0);

//if (m.is_ring_atom (zatom))
//  _do_not_process_rings_containing (m, zatom, tkfta);

  return;
}

int
Toggle_Kekule_Form::process (Molecule & m, 
                             const Set_of_Atoms & embedding,
                             int & changed)
{
  int nb = _bond.number_elements ();
  if (0 == nb)
  {
    cerr << "Toggle_Kekule_Form::process: no bonds specified\n";
    return 0;
  }

  if (! ok_embedding (embedding))
  {
    cerr << "Toggle_Kekule_Form::process: invalid embedding " << embedding << endl;
    return 0;
  }

  if (_all_bonds_correct (m, embedding))
  {
    changed = 0;
    return 1;
  }

// Looks like we are going to have to change the molecule

  m.compute_aromaticity ();
  if (! _all_bonds_aromatic (m, embedding))
  {
    cerr << "Toggle_Kekule_Form::process:not all bonds aromatic '" << m.name () << "'\n";
    return 0;
  }

// All the methods need a means of determining whether or not the
// bond between any two atoms can be changed.

  Toggle_Kekule_Form_Temporary_Arrays tkfta (m);

  _do_chemistry (m, tkfta);

#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
  cerr << "Aft4er _do_chemistry\n";
  for (int i = 0; i < m.natoms(); i++)
  {
    for (int j = i + 1; j < m.natoms(); j++)
    {
      if (! m.are_bonded(i, j))
        continue;

      if (! m.in_same_ring(i, j))
        continue;

//    if (can_change_bond[i * matoms + j])
        cerr << "Can change bond between atoms " << i << " and " << j << " value " << tkfta.can_change_bond()[i*m.natoms()+j] << endl;
    }
  }
#endif

  int rc = _process (m, embedding, tkfta);

  if (rc)
  {
    changed = 1;
    m.invalidate_smiles ();
  }

  return rc;
}

/*
  Note that outside_ring may be a ring atom - just not part of the
  set of atoms being processed
*/

/*static boolean
multiple_bond_outside_ring (Molecule & m,
                            atom_number_t i, 
                            const Atom * a,
                            int id,
                            const int * process_these)
{
  int acon = a->ncon ();
  for (int j = 0; j < acon; j++)
  {
    const Bond * b = a->item (j);
    if (b->is_single_bond ())
      continue;

    atom_number_t k = b->other (i);

    if (id != process_these[k])
      return true;
  }

  return false;
}*/

static boolean
multiple_bond_outside_aromatic_ring (Molecule & m,
                                     atom_number_t i, 
                                     const Atom * a)
{
  int acon = a->ncon ();
  for (int j = 0; j < acon; j++)
  {
    const Bond * b = a->item (j);
    if (b->is_single_bond ())
      continue;

//  Got a double bond

    atom_number_t k = b->other (i);

    if (m.is_non_ring_atom (k))    // to a non ring atom
      return true;

    if (m.is_aromatic (k))    // either in the ring our part of fused system
      continue;

    return true;     // a non-aromatic ring atom
  }

  return false;
}

//#define DEBUG_GROW_RING_SYSTEM

static int
grow_ring_system (Molecule & m,
                  const Ring * r,
                  int id,
                  Toggle_Kekule_Form_Temporary_Arrays & tkfta,
                  int * ring_already_done)
{
  int rc = 1;

  r->set_vector (tkfta.process_these (), id);

  int frn = r->fused_ring_neighbours ();

#ifdef DEBUG_GROW_RING_SYSTEM
  cerr << "grow_ring_system: rid " << id  << ", ring number " << r->ring_number () << ", " << frn << " fused neighbours\n";
#endif

  for (int j = 0; j < frn; j++)
  {
    const Ring * rj = r->fused_neighbour (j);

    int rjn = rj->ring_number ();

#ifdef DEBUG_GROW_RING_SYSTEM
    cerr << "ring " << r->ring_number () << " joined to ring " << rjn << endl;
    assert (rjn >= 0 && rjn < m.nrings ());
#endif

    if (ring_already_done[rjn])
      continue;

    if (! rj->is_aromatic ())
      continue;

    if (! tkfta.ring_can_toggle (rjn))
      continue;

    ring_already_done[rjn] = 1;
    rc += grow_ring_system (m, rj, id, tkfta, ring_already_done);
  }

#ifdef DEBUG_GROW_RING_SYSTEM
  cerr << "grow_ring_system returns " << rc << endl;
#endif

  return rc;
}

int
Toggle_Kekule_Form::_process (Molecule & m,
                              const Set_of_Atoms & embedding,
                              Toggle_Kekule_Form_Temporary_Arrays & tkfta)
{
  int matoms = m.natoms ();

  int * can_change_bond = tkfta.can_change_bond ();

//#define DEBUG_DO_TOGGLE_KEKULE_FORM
#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
  cerr << "Start of _process\n";
  for (int i = 0; i < matoms; i++)
  {
    for (int j = i + 1; j < matoms; j++)
    {
      if (! m.are_bonded(i, j))
        continue;

      if (! m.in_same_ring(i, j))
        continue;

//    if (can_change_bond[i * matoms + j])
        cerr << "Can change bond between atoms " << i << " and " << j << " value " << can_change_bond[i*matoms+j] << endl;
    }
  }
#endif

  int nb = _bond.number_elements ();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond[i];

    if (0 == can_change_bond[b->a1() * matoms + b->a2()])   // we cannot do anything
      return 0;

    atom_number_t a1 = embedding[b->a1 ()];
    atom_number_t a2 = embedding[b->a2 ()];

    if (_correct[i])
      can_change_bond[a1 * matoms + a2] = can_change_bond[a2 * matoms + a1] = 0;
  }

  int * process_these = tkfta.process_these ();

// figure out which rings are involved

  int nr = m.nrings ();

  if (1 == nr)
  {
    const Ring * r = m.ringi (0);

    r->set_vector (process_these, 0);

    return _process_single_ring (m, embedding, 0, tkfta);
  }

// We may have fused rings
// We want to set up the process_these array before changing the molecule.

  resizable_array<int> isolated_rings, fused_rings;

  isolated_rings.resize (nr);
  fused_rings.resize (nr);

  int * ring_already_done = new_int (nr); iw_auto_array<int> free_ring_already_done (ring_already_done);

#ifdef DEBUG_GROW_RING_SYSTEM
  for (int i = 0; i < nr; i++)
  {
    cerr << "Ring " << i << " can toggle? " << tkfta.ring_can_toggle (i) << endl;
  }
#endif

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])      // part of a fused system done before
      continue;

    if (! tkfta.ring_can_toggle (i))
      continue;

    const Ring * r = m.ringi (i);

    if (r->is_non_aromatic ())     // very important!
      continue;

    if (! will_change_ring (r, embedding))    // none of the bonds we insist upon are in the ring
      continue;

    int rid = i + 1;     // make sure not zero (since +0 == -0)

    if (! r->is_fused ())
    {
      r->set_vector (process_these, rid);

      isolated_rings.add (rid);
    }
    else if (1 == grow_ring_system (m, r, rid, tkfta, ring_already_done))
    {
      isolated_rings.add (rid);
    }
    else 
      fused_rings.add (rid);

    ring_already_done[i] = 1;    // probably not necessary
  }

#ifdef DEBUG_GROW_RING_SYSTEM
  cerr << isolated_rings.number_elements () << " isolated and " << fused_rings.number_elements () << " fused rings\n";
#endif

#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
  cerr << "Processing " << isolated_rings.number_elements () << " isolated rings\n";
#endif

  int rc = 1;
  for (int i = 0; i < isolated_rings.number_elements (); i++)
  {
    int rid = isolated_rings[i];
    if (! _process_single_ring (m, embedding, rid, tkfta))
    {
      rc = 0;
      if (_display_error_messages)
        cerr << "Toggle_Kekule_Form::_process: cannot change isolated rid " << rid << ", molecule '" << m.name () << "'\n";
    }
  }

#ifdef DEBUG_DO_TOGGLE_KEKULE_FORM
  cerr << "Processing " << fused_rings.number_elements () << " fused ring sets\n";
#endif

  for (int i = 0; i < fused_rings.number_elements (); i++)
  {
    int rid = fused_rings[i];
    if (! _process_ring_system (m, embedding, rid, tkfta))
    {
      rc = 0;
      if (_display_error_messages)
        cerr << "Toggle_Kekule_Form::_process: cannot change fused rid " << rid << ", molecule '" << m.name () << "'\n";
    }
  }

  return rc;
}

void
Toggle_Kekule_Form::_do_chemistry_aromatic_ring (Molecule & m,
                                                 const Ring & r,
                                                 Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  Atom * const * atom = tkfta.atom ();

  int ring_size = r.number_elements ();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

    Atom * a = atom[j];
    
    atomic_number_t z = a->atomic_number ();

    int acon = a->ncon ();

    if (6 == z && 2 == acon)     // these will vary
      continue;

    int nb = a->nbonds ();

    if (7 == z && 2 == acon && 3 == nb)    // like pyridine
      continue;

    if (7 == z && 2 == acon && 2 == nb)    // like in N1C=CC=C1, pyrrole
    {
      if (! _allow_pyrrole_to_change)
        _no_changes_to_atom (m, j, tkfta);
      continue;
    }

    if (7 == z && 3 == acon && 3 == nb)    // methyl-pyrrole
    {
      _no_changes_to_atom (m, j, tkfta);
      continue;
    }

    if (8 == z || 16 == z)
    {
      _no_changes_to_atom (m, j, tkfta);
      continue;
    }

    if (6 == z && 3 == acon && 4 == nb && multiple_bond_outside_aromatic_ring (m, j, a))
    {
      _no_changes_to_atom (m, j, tkfta);
      continue;
    }
  }

  return;
}

void
Toggle_Kekule_Form::_do_chemistry (Molecule & m,
                                   Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  int nr = m.nrings ();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi (i);

    if (ri->is_aromatic ())
      _do_chemistry_aromatic_ring (m, *ri, tkfta);
    else
      tkfta.set_ring_can_toggle (i, 0);
  }

  return;
}

//#define DEBUG_NEXT_ATOM_IN_RING

static atom_number_t
next_atom_in_ring (const Atom * a,
                   int rid,
                    const int * process_these,
                    atom_number_t aprev,
                    atom_number_t astart2)
{
  int acon = a->ncon ();

#ifdef DEBUG_NEXT_ATOM_IN_RING
  cerr << "Looking for next atom in ring from " << astart2 << ", id = " << rid << endl;
#endif

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (astart2, i);

    if (j == aprev)    // need to move around the ring
      continue;

    if (rid == process_these[j])
      return j;
  }

  assert (NULL == "Did not find next atom in ring for Kekule Toggle");

  return INVALID_ATOM_NUMBER;
}

void
Toggle_Kekule_Form::_set_all_bonds_to_single (Molecule & m,
                                              int id,
                                              Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  const int * process_these = tkfta.process_these ();
  const int * can_change_bond = tkfta.can_change_bond ();
  Atom * const * atom = tkfta.atom ();

  int matoms = m.natoms ();

  for (int i = 0; i < matoms; i++)
  {
    if (id != process_these[i])
      continue;

    Atom * a = atom[i];

    int acon = a->ncon ();

    if (acon == a->nbonds ())      // atom already has single bonds.
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item (j);
      if (b->is_single_bond ())
        continue;

      atom_number_t k = b->other (i);

      if (k < i)
        continue;

      if (id != process_these[k])
        continue;

//    cerr << "Can we change bond between " << i << " and " << k << " " << can_change_bond[i*matoms+k] << endl;
      if (0 == can_change_bond[i * matoms + k])
        continue;

      m.set_bond_type_between_atoms (i, k, SINGLE_BOND);
    }
  }
}

int
Toggle_Kekule_Form::_set_our_bonds (Molecule & m,
                                    const Set_of_Atoms & embedding,
                                    int id,
                                    Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  const int * process_these = tkfta.process_these ();
  int * can_change_bond = tkfta.can_change_bond ();
  int * has_double_bond = tkfta.has_double_bond ();

  int matoms = m.natoms ();

  for (int i = 0; i < _bond.number_elements (); i++)
  {
    if (_correct[i])
      continue;

    const Bond * b = _bond[i];

    atom_number_t a1 = embedding[b->a1 ()];
    atom_number_t a2 = embedding[b->a2 ()];

    if (id != process_these[a1] || id != process_these[a2])
      continue;

    can_change_bond[a1 * matoms + a2] = can_change_bond[a2 * matoms + a1] = 0;

    if (b->is_double_bond ())
    {
      m.set_bond_type_between_atoms (a1, a2, DOUBLE_BOND);
      has_double_bond[a1] = 2;
      has_double_bond[a2] = 2;
    }
  }

  return 1;
}

int
Toggle_Kekule_Form::_process_single_ring (Molecule & m,
                                          const Set_of_Atoms & embedding,
                                          int rid,
                                          Toggle_Kekule_Form_Temporary_Arrays & tkfta)
{
  const int * process_these = tkfta.process_these();
  Atom * const * atom = tkfta.atom ();

  resizable_array_p<Bond> existing_bonds;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (rid != process_these[i])
      continue;

    const Atom * ai = atom[i];

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item(j);

      atom_number_t k = b->other(i);

      if (k < i)
        continue;

      if (b->is_single_bond())
        existing_bonds.add(new Bond (i, k, SINGLE_BOND));
      else
        existing_bonds.add(new Bond (i, k, DOUBLE_BOND));
    }
  }

  if (_process_single_ring2(m, embedding, rid, tkfta))
    return 1;

  for (int i = existing_bonds.number_elements() - 1; i >= 0; i--)
  {
    const Bond * b = existing_bonds[i];

    m.set_bond_type_between_atoms(b->a1(), b->a2(), b->btype());
  }

  return 0;
}

int
Toggle_Kekule_Form::_process_single_ring2 (Molecule & m,
                                           const Set_of_Atoms & embedding,
                                           int rid,
                                           Toggle_Kekule_Form_Temporary_Arrays & tkfta)
{
  atom_number_t astart1 = embedding[_bond[0]->a1 ()];
  atom_number_t astart2 = embedding[_bond[0]->a2 ()];

  if (! tkfta.atom_can_change(astart1) || ! tkfta.atom_can_change(astart2))
    return 0;

  _set_all_bonds_to_single (m, rid, tkfta);

  _set_our_bonds (m, embedding, rid, tkfta);

  int matoms = m.natoms ();

  assert (rid >= 0);

  int * can_change_bond = tkfta.can_change_bond ();
  int * process_these = tkfta.process_these ();
  Atom * const * atom = tkfta.atom ();

// We need to keep track of our starting point (astart1) and the atom being
// processed (astart2)

  bond_type_t previous_bond = _bond[0]->btype ();

//#define DEBUG_TOGGLE_KEKULE_PROCESS_SINGLE_RING
#ifdef DEBUG_TOGGLE_KEKULE_PROCESS_SINGLE_RING
  cerr << "Starting smiles '" << m.smiles () << "'\n";
  cerr << "Start with atoms " << astart1 << " and " << astart2 << endl;
  cerr << "Initial bond is " << previous_bond << endl;
#endif

  atom_number_t aprev = astart1;
  while (1)
  {
    atom_number_t anext = next_atom_in_ring (atom[astart2], rid, process_these, aprev, astart2);

    if (DOUBLE_BOND == previous_bond)    // cannot have consecutive double bonds
      previous_bond = SINGLE_BOND;
    else if (! can_change_bond[astart2 * matoms + anext])    //
      previous_bond = BOND_TYPE_ONLY (m.bond_between_atoms (astart2, anext)->btype ());
    else if (0 == m.hcount(anext))   // no room for a double bond
      previous_bond = BOND_TYPE_ONLY (m.bond_between_atoms (astart2, anext)->btype ());
    else
    {
      m.set_bond_type_between_atoms (astart2, anext, DOUBLE_BOND);
      previous_bond = DOUBLE_BOND;
    }

#ifdef DEBUG_TOGGLE_KEKULE_PROCESS_SINGLE_RING
    cerr << "Atom pair is " << astart2 << " to " << anext << " set to " << previous_bond << " smiles '" << m.smiles () << "'\n";
#endif

    if (astart1 == anext)     // we have looped back to the beginning
      break;

    aprev = astart2;
    astart2 = anext;
  }

#ifdef DEBUG_TOGGLE_KEKULE_PROCESS_SINGLE_RING
  cerr << "Final smiles '" << m.smiles () << "'\n";
#endif

  return _all_atoms_aromatic (m, rid, _display_error_messages, tkfta);
}

int
Toggle_Kekule_Form::_get_ring_system_atoms (resizable_array<int> & atoms_to_process,
                       int rid,
                       atom_number_t zatom,
                       Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  const Atom * a = tkfta.atom ()[zatom];
  int * process_these = tkfta.process_these ();

  int rc = 0;

  int acon = a->ncon ();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (zatom, i);

    if (rid != process_these[j])
      continue;

    atoms_to_process.add (j);
    process_these[j] = -rid;

    rc++;

//  cerr << "From atom " << zatom << " go to atom " << j << endl;
    rc += _get_ring_system_atoms (atoms_to_process, rid, j, tkfta);
  } 

  return rc;
}

int
Toggle_Kekule_Form::_process_ring_system (Molecule & m,
                                          const Set_of_Atoms & embedding,
                                          int rid,
                                          Toggle_Kekule_Form_Temporary_Arrays & tkfta)
{
//_do_chemistry (m, rid, tkfta);

  _set_all_bonds_to_single (m, rid, tkfta);

  _set_our_bonds (m, embedding, rid, tkfta);

  int * process_these = tkfta.process_these ();

// Find a starting point with a known bond

  atom_number_t astart1 = INVALID_ATOM_NUMBER;
  atom_number_t astart2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _bond.number_elements (); i++)
  {
    atom_number_t a1 = embedding[_bond[i]->a1 ()];
    if (rid != process_these[a1])
      continue;

    if (! tkfta.atom_can_change(a1))
      return 0;

    atom_number_t a2 = embedding[_bond[i]->a2 ()];
    if (rid != process_these[a2])
      continue;

    if (! tkfta.atom_can_change(a2))
      return 0;

    astart1 = a1;
    astart2 = a2;
    break;
  }

  assert (INVALID_ATOM_NUMBER != astart1);

#ifdef DEBUG_PROCESS_RING_SYSTEM
  cerr << "Starting with atoms " << astart1 << " and " << astart2 << endl;
#endif

  resizable_array<atom_number_t> atoms_to_process;

  int matoms = m.natoms ();

  atoms_to_process.resize (matoms);

  atoms_to_process.add (astart1);
  process_these[astart1] = -rid;
  atoms_to_process.add (astart2);
  process_these[astart2] = -rid;

  _get_ring_system_atoms (atoms_to_process, rid, astart1, tkfta);

  for (int i = 0; i < matoms; i++)   // turn back to positive numbers
  {
    if (-rid == process_these[i])
      process_these[i] = rid;
  }

  int rc = _process_ring_system (m, atoms_to_process, rid, 0, tkfta);

  return rc;
}

int
Toggle_Kekule_Form::_process_ring_system (Molecule & m,
                                          resizable_array<int> & atoms_to_process,
                                          int rid,
                                          int zitem,
                                          Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  atom_number_t curr = atoms_to_process[zitem];

  int * process_these = tkfta.process_these ();
  int * has_double_bond = tkfta.has_double_bond ();
  Atom * const * atom = tkfta.atom ();

#ifdef DEBUG_PROCESS_RING_SYSTEM
  cerr << "Processing ring system " << rid << endl;
  for (int i = 0; i < m.natoms (); i++)
  {
    if (rid == process_these[i])
      cerr << "Processing atom " << i << " " << m.smarts_equivalent_for_atom (i) << endl;
  }
#endif

  process_these[curr] = -rid;     // so we don't double back on ourselves

  const Atom * acurr = atom[curr];

  int acon = acurr->ncon ();

// First decide whether or not this atom can get a double bond within the ring

  int tochange;

  if (! tkfta.atom_can_change (curr))
    tochange = 0;
  else if (acurr->nbonds () > acon)     // what about *=N(=O)-*
    tochange = 0;
  else
    tochange = 1;

#ifdef OLD_CODE
// Seems as if all these checks can be replace by the unsaturation test above
  else if (has_double_bond[curr])     // already has a double bond - what about *=N(=O)-* groups?
    tochange = 0;
  else if (6 == z && acurr->nbonds () > acon)      // carbon's can have only 1 double bond
    tochange = 0;
  else if (7 == z && 2 == acon && acurr->nbonds () > acon)    // like pyridine
    tochange = 0;
  else if (7 == z && 3 == acon && 5 == acurr->nbonds ())
    tochange = 0;
#endif

#ifdef DEBUG_PROCESS_RING_SYSTEM
  cerr << "Atom " << curr << " '" << m.smarts_equivalent_for_atom (curr) << " tochange? " << tochange << endl;
#endif

  if (! tochange)
  {
    int rc;
    if (zitem == atoms_to_process.number_elements () - 1)    // no more to test
      rc = _all_atoms_aromatic (m, -rid, 0, tkfta);
    else
      rc = _process_ring_system (m, atoms_to_process, rid, zitem+1, tkfta);

    process_these[curr] = rid;

    return rc;
  }

  int matoms = m.natoms ();

// Somewhere there must be a connection where a double bond makes sense

  const int * can_change_bond = tkfta.can_change_bond ();

  for (int i = 0; i < acon; i++)
  {
    const Bond * existing_bond = acurr->item (i);

    if (existing_bond->is_double_bond ())     // we are hoping to place a double bond
      continue;

    atom_number_t j = existing_bond->other (curr);

    if (rid != process_these[j])
      continue;

    if (0 == can_change_bond[curr * matoms + j])
      continue;

    if (has_double_bond[j])
      continue;

    m.set_bond_type_between_atoms (curr, j, DOUBLE_BOND);
    has_double_bond[curr] = 1;
    has_double_bond[j] = 1;

#ifdef DEBUG_PROCESS_RING_SYSTEM
    cerr << "Set double bond between atoms " << curr << " and " << j << " '" << m.smiles () << "'\n";
#endif

    int rc;

    if (zitem == atoms_to_process.number_elements () - 1)
      rc = _all_atoms_aromatic (m, -rid, 0, tkfta);
    else
      rc = _process_ring_system (m, atoms_to_process, rid, zitem + 1, tkfta);

    if (0 == rc)    // placing the double bond here didn't work, back to single bond
    {
#ifdef DEBUG_PROCESS_RING_SYSTEM
      cerr << "Reset bond between " << curr << " and " << j << " back to single\n";
#endif
      m.set_bond_type_between_atoms (curr, j, SINGLE_BOND);
      has_double_bond[curr] = 0;
      has_double_bond[j] = 0;
    }

    if (rc)
      return rc;
  }

  process_these[curr] = rid;

  return 0;     // no kekule form found
}

/*
  We have finished a ring or ring-system. All the atoms in that
  system should now be aromatic.
*/

int
Toggle_Kekule_Form::_all_atoms_aromatic (Molecule & m,
                                         int id,
                                         int display_error_message,
                                         Toggle_Kekule_Form_Temporary_Arrays & tkfta) const
{
  const int * process_these = tkfta.process_these ();
  aromaticity_type_t * arom = tkfta.arom ();

  m.invalidate_smiles ();

  m.aromaticity (arom);

  int rc = 1;

  int matoms = m.natoms ();

  for (int i = 0; i < matoms; i++)
  {
    if (id != process_these[i])
      continue;

    if (! arom[i])
    {
      if (display_error_message)
        cerr << "Toggle_Kekule_Form::_process_single_ring: atom " << i << " not aromatic (" << m.atomic_symbol (i) << ' ' << m.ncon (i) << " connections)\n";
      rc = 0;
    }
  }

  if (0 == rc && display_error_message)
    cerr << "Smiles is '" << m.smiles () << "'\n";

  return rc;
}

int
Toggle_Kekule_Form::add_bond (int a1,
                              int a2,
                              bond_type_t bt)
{
  assert (a1 >= 0 && a2 >= 0 && a1 != a2);

  Bond * b;
  if (IS_SINGLE_BOND (bt))
    b = new Bond (a1, a2, SINGLE_BOND);
  else if (IS_DOUBLE_BOND (bt))
    b = new Bond (a1, a2, DOUBLE_BOND);
  else if (IS_TRIPLE_BOND (bt))
    b = new Bond (a1, a2, TRIPLE_BOND);
  else
  {
    cerr << "Toggle_Kekule_Form::add_bond:unrecognised bond type " << bt << endl;
    return 0;
  }

  return add_bond (b);
}

int
Toggle_Kekule_Form::add_bond (Bond * b)
{
  _bond.add (b);

  if (NULL != _correct)
    delete [] _correct;

  _correct = new int[_bond.number_elements ()];

//cerr << "_correct allocated to " << _bond.number_elements () << " items\n";

  return 1;
}

int
Toggle_Kekule_Form::add_bond_from_msi_attribute (const msi_attribute & msi)
{
  if (3 != msi.number_int_values ())
  {
    cerr << "Toggle_Kekule_Form::add_bond_from_msi_attribute:invalid specification '" << msi << "'\n";
    return 0;
  }

  int a1 = msi.int_multi_value (0);
  int a2 = msi.int_multi_value (1);

  if (a1 < 0 || a2 < 0 || a1 == a2)
  {
    cerr << "Toggle_Kekule_Form::add_bond_from_msi_attribute:invalid atom specification '" << msi << "'\n";
    return 0;
  }

  int b = msi.int_multi_value (2);
  bond_type_t bt;

  if (1 == b)
    bt = SINGLE_BOND;
  else if (2 == b)
    bt = DOUBLE_BOND;
  else
  {
    cerr << "Toggle_Kekule_Form::add_bond_from_msi_attribute:unrecognised bond type '" << msi << "'\n";
    return 0;
  }

  return add_bond (new Bond (a1, a2, bt));
}

Toggle_Kekule_Form_Temporary_Arrays::Toggle_Kekule_Form_Temporary_Arrays (Molecule & m)
{
  int matoms = m.natoms ();

  _can_change_bond = new_int (matoms * matoms, 1);

  _process_these = new_int (matoms, numeric_limits<int>::min());

  _arom = new aromaticity_type_t[matoms];

  m.aromaticity (_arom);

  _atom = new Atom *[matoms];

  m.atoms ( (const Atom **) _atom);

  _has_double_bond = new_int (matoms, 0);

  _ring_can_vary = new_int (m.nrings (), 1);

  _atom_can_change = new_int (matoms, 1);

  return;
}

Toggle_Kekule_Form_Temporary_Arrays::~Toggle_Kekule_Form_Temporary_Arrays ()
{
  delete _can_change_bond;
  delete _process_these;
  delete _arom;
  delete _atom;
  delete _has_double_bond;
  delete _ring_can_vary;
  delete _atom_can_change;
  
  return;
}

static int
btype_back_to_mdl_form (bond_type_t b)
{
  if (SINGLE_BOND == b)
    return 1;

  if (DOUBLE_BOND == b)
    return 2;

  if (TRIPLE_BOND == b)
    return 3;

  cerr << "btype_back_to_mdl_form: warning, unimplemented feature " << b << endl;

  return 5;
}

int
Toggle_Kekule_Form::write_msi (ostream & os,
                               const IWString & ind,
                               const char * attribute_name) const
{
  for (int i = 0; i < _bond.number_elements (); i++)
  {
    const Bond * b = _bond[i];

    os << ind << "(A I " << attribute_name << " (" << b->a1 () << ' ' << b->a2 () << ' ' << btype_back_to_mdl_form (b->btype ()) << "))\n";
  }

  return os.good ();
}

/*
  this is kind of bizzare. All we do is set up the conditions where the regular process() member
  function can be called...
*/

int
Toggle_Kekule_Form::process (Molecule & m,
                             atom_number_t a1,
                             atom_number_t a2,
                             bond_type_t bt,
                             int & changed)
{
  assert (m.are_bonded(a1, a2));

  const Bond * b = m.bond_between_atoms(a1, a2);

  assert (NULL != b);

  changed = 0;

  if (SINGLE_BOND == bt && b->is_single_bond())
    return 1;

  if (DOUBLE_BOND == bt && b->is_double_bond())
    return 1;

  Set_of_Atoms s;
  s.add(a1);
  s.add(a2);

  Bond * dummy = new Bond(0, 1, bt);

  _bond.resize_keep_storage(0);

  _bond.add(dummy);

  if (NULL != _correct)
    delete [] _correct;

  _correct = new_int (1); 

  int rc = process(m, s, changed);

  delete [] _correct;
  _correct = NULL;

  _bond.resize_keep_storage(0);

  return rc;
}
