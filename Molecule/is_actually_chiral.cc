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
#include <iostream>

using namespace std;

#include "misc.h"

#include "molecule.h"
#include "path_scoring.h"
#include "chiral_centre.h"
#include "is_actually_chiral.h"

static int max_iterations = 0;

void
set_max_iterations (int m)
{
  assert (m > 0);

  max_iterations = m;
}

/*
  To determine if an atom is chiral or not, we need to perform path tracing
  from that atom.
*/

static int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom,
                    resizable_array_p<Path_Scoring> & ps,
                    int * claimed,
                    Atom * const * atom)
{
  const Atom * a = atom[zatom];

  int acon = a->ncon ();

  if (ps.number_elements ())
    ps.resize_keep_storage (0);

  ps.resize (acon);

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (zatom, i);

    Path_Scoring * p = new Path_Scoring;

    const Atom * aj = atom[j];
    p->initialise (j, aj);

    claimed[j] = 1;

    ps.add (p);
  }

  int stopped;
  if (resolved (ps, stopped))
    return 1;

  int iterations = 0;

  while (1)
  {
    for (int i = 0; i < acon; i++)
    {
      if (ps[i]->active ())
        ps[i]->advance (atom, claimed);
    }

    int stopped;
    if (resolved (ps, stopped))
      return 1;

    if (stopped)      // not resolved, but cannot go any further
      return 0;

    int number_active = 0;
    for (int i = 0; i < acon; i++)
    {
      if (! ps[i]->active ())
        continue;

      ps[i]->update_claimed (claimed);
      number_active++;
    }

    if (number_active < 2)
      return 0;

    iterations++;

    if (max_iterations > 0 && iterations >= max_iterations)
    {
      cerr << "Not resolved by " << max_iterations << " iterations\n";
      return 0;
    }
  }

  return 1;
}

/*
  The query for an asymmetric carbon atom will hit things like the
  carbon in t-butyl. We need to examine the neighbours to make sure
  that this atom actually is an asymmetric centre

  Thought about putting in a more aggressive check on the number of
  connections, but too dangerous. Even 2 == ncon is problematic because
  you could have an atom with a lone-pair and an implicit Hydrogen
*/

int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom)
{
  resizable_array_p<Path_Scoring> ps;

  return is_actually_chiral (m, zatom, ps);
}

int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom,
                    resizable_array_p<Path_Scoring> & ps)
{
  if (1 == m.ncon (zatom))
    return 0;

  if (m.ncon (zatom) > 4)
    return 0;

  if (m.hcount (zatom) > 1)    // what if isotopic Hydrogen???
    return 0;

  int lp;
  if (m.lone_pair_count (zatom, lp) && lp > 1)
    return 0;

  m.compute_aromaticity_if_needed ();    // so bonds get aromatic character

  int matoms = m.natoms ();

  int * claimed = new_int (matoms);

  claimed[zatom] = 1;

  Atom * const * atoms = new Atom *[matoms];

  m.atoms ( (const Atom **) atoms);

  int rc = is_actually_chiral (m, zatom, ps, claimed, atoms);

  delete claimed;
  delete atoms;

  return rc;
}

int
do_remove_invalid_chiral_centres (Molecule & m)
{
  int nc = m.chiral_centres ();
  if (0 == nc)
    return 0;

// Removing a chiral centre while we are scanning the set would mess things up,
// so we make a list of the atoms with invalid chiral centres and remove them later

  Set_of_Atoms centres_to_be_removed;

  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number (i);

    atom_number_t a = c->a ();

    if (! is_actually_chiral (m, a))
      centres_to_be_removed.add (a);
  }

  if (centres_to_be_removed.number_elements ())
  {
    for (int i = 0; i < centres_to_be_removed.number_elements (); i++)
    {
      m.remove_chiral_centre_at_atom (centres_to_be_removed[i]);
    }
  }

  return nc;
}
