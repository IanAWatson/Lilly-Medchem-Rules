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
  Implementations of atom object functions.
*/

#include <iostream>
#include <iomanip>
#if (__GNUC__ == 2)
#include <function.h>
#else 
#include <functional>
#endif

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

using namespace std;

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "atom.h"
#include "misc2.h"
#include "iwrandom.h"

/*
  This is a count of the number of atoms created - the number
  freed.
*/

static int current_atoms = 0;

int
how_many_atoms()
{
  return current_atoms;
}

static int copy_implicit_hydrogen_count_in_atom_copy_constructor = 1;

void
set_copy_implicit_hydrogen_count_in_atom_copy_constructor (int s)
{
  copy_implicit_hydrogen_count_in_atom_copy_constructor = s;
}

static formal_charge_t _min_reasonble_atomic_formal_charge_value = -5;
static formal_charge_t _max_reasonble_atomic_formal_charge_value =  5;

int
set_reasonable_formal_charge_range (formal_charge_t cmin, formal_charge_t cmax)
{
  _min_reasonble_atomic_formal_charge_value = cmin;
  _max_reasonble_atomic_formal_charge_value = cmax;

  return 1;
}

int
reasonable_formal_charge_value (formal_charge_t c)
{
  if (c < _min_reasonble_atomic_formal_charge_value)
    return 0;

  if (c > _max_reasonble_atomic_formal_charge_value)
    return 0;

  return 1;
}

static int alternate_valences_give_hcount = 1;

void 
set_alternate_valences_give_hcount(const int s)
{
  alternate_valences_give_hcount = s;

  return;
}


static int file_scope_four_connected_neutral_nitrogen_has_h = 0;

void
set_four_connected_neutral_nitrogen_has_h(const int s)
{
  file_scope_four_connected_neutral_nitrogen_has_h = s;
}

void
Atom::_default_values (const Element * zelement)
{
  assert (OK_ELEMENT (zelement));

//current_atoms++;

  _element = zelement;

  _formal_charge = 0;

  _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;

  _implicit_hydrogens_known = 0;

  _nrings = ATOM_PROPERTY_UNKNOWN;

  _nbonds = ATOM_PROPERTY_UNKNOWN;

  _permanent_aromatic = zelement->permanent_aromatic();

  _isotope = 0;

  _user_specified_void_ptr = NULL;

  resize(3);   // waste some space for efficiency. In general, this is good

  return;
}

Atom::Atom (const Element * zelement)
{
  _default_values(zelement);

  return;
}

Atom::Atom (const char * asymbol)
{
  assert (NULL != asymbol);

  const Element * e = get_element_from_symbol(asymbol, static_cast<int>(::strlen(asymbol)), _isotope);

  if (NULL == e)
  {
    if (auto_create_new_elements())
      e = new Element(asymbol);

    if (NULL == e)
    {
      cerr << "Atom::Atom unrecognised atomic symbol '" << asymbol << "\n";
      iwabort();
    }
  }

  _default_values(e);

  return;
}

Atom::Atom (atomic_number_t zz)
{
  const Element * e = get_element_from_atomic_number(zz);
  if (NULL == e)
  {
    cerr << "Atom::Atom: cannot get element for z = " << zz << endl;
    iwabort();
  }

  _default_values(e);

  return;
}

Atom::~Atom()
{
  assert(ok());

//_x = _y = _z = coord_t(0.0);

//current_atoms--;

  return;
}

/*
  Make a copy of another atom. Connections not copied.
*/

void
Atom::_constructor_copy_atom_attributes (const Atom & other_atom)
{
  assert(other_atom.ok());

  _default_values(other_atom._element);

  _x = other_atom._x;
  _y = other_atom._y;
  _z = other_atom._z;

  _implicit_hydrogens_known = other_atom._implicit_hydrogens_known;

  if (_implicit_hydrogens_known)
    _implicit_hydrogens = other_atom._implicit_hydrogens;
  else if (copy_implicit_hydrogen_count_in_atom_copy_constructor)
    _implicit_hydrogens = other_atom._implicit_hydrogens;
  else
    _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;

  _isotope = other_atom._isotope;

  _formal_charge = other_atom._formal_charge;

  _user_specified_void_ptr = other_atom._user_specified_void_ptr;

  resize(other_atom._elements_allocated);

  return;
}

Atom::Atom (const Atom & other_atom)
{
  _constructor_copy_atom_attributes(other_atom);

  return;
}

Atom::Atom (const Atom * other_atom)
{
  _constructor_copy_atom_attributes(*other_atom);

  return;
}

int
Atom::debug_print (ostream & os) const
{
  assert (os.good());
  
  os << "Atom details " << this << endl;

  if (! OK_ELEMENT(_element))
    os << "Our element pointer seems to be messed up\n";
  else
    _element->debug_print(os);

  os << "Ncon = " << number_elements() << endl;
  for (int i = 0; i < number_elements(); i++)
  {
    item (i)->debug_print (os);
    os << endl;
  }

  if (ATOM_PROPERTY_UNKNOWN != _implicit_hydrogens)
  {
    os << "hcount " << _implicit_hydrogens;
    if (_implicit_hydrogens_known)
      os << '*';
    os << endl;
  }

  os.setf(std::ios::showpoint);
  os << "Coordinates (" << _x << "," << _y << "," << _z << ")\n";

  return 1;
}

/*
  If we are getting an extra connection, we also need set_modified
*/

int
Atom::add (Bond * b)
{
  resizable_array<Bond *>::add(b);

  set_modified();

  return 1;
}









int
Atom::nbonds()
{
  assert (ok());   

  if (ATOM_PROPERTY_UNKNOWN == _nbonds)
  {
    _nbonds = 0;
    for (int i = 0; i < _number_elements; i++)
    {
      _nbonds += _things[i]->number_of_bonds();
    }
  }

  return _nbonds;
}

/*
  A const version
*/

int
Atom::nbonds() const
{
  assert (ok());   

  if (ATOM_PROPERTY_UNKNOWN != _nbonds)
    return _nbonds;

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->number_of_bonds();
  }

  return rc;
}

int
Atom::fully_saturated() const
{
  equal_to<int> c;
  return _common_saturation(c);
}

int
Atom::unsaturated() const
{
  less<int> l;
  return _common_saturation(l);
}

template <typename T>
int
Atom::_common_saturation(const T & comparitor) const
{
  if (ATOM_PROPERTY_UNKNOWN == _nbonds)
  {
    int nb = 0;
    for (int i = 0; i < _number_elements; i++)
    {
      nb += _things[i]->number_of_bonds();
    }

    return comparitor(nb, _number_elements);
  }

  return comparitor(_nbonds, _number_elements);
}

template int Atom::_common_saturation(const equal_to<int> &) const;
template int Atom::_common_saturation(const less<int> &) const;

int
Atom::recompute_nbonds()
{
  _nbonds = ATOM_PROPERTY_UNKNOWN;

  return nbonds();
}

const IWString &
Atom::atomic_symbol() const
{
  assert (ok());

  return _element->symbol();
}

int
Atom::atomic_number() const
{
  assert (ok());

  return _element->atomic_number();
}

atomic_mass_t
Atom::atomic_weight() const
{
  assert (ok());

  return _element->atomic_mass();
}

int
Atom::exact_mass (exact_mass_t & zresult) const
{
  assert (ok());

  zresult = _element->exact_mass();

  if (0.0 == zresult)
    return 0;

  return 1;
}

const Element *
Atom::element() const
{
  assert (ok());

  return _element;
}

const Element &
Atom::elementq() const
{
  assert (ok());

  return *_element;
}

void
Atom::set_element (const Element * new_element)
{
  assert (ok());
  assert (OK_ELEMENT (new_element));

  if (_element && (_element->atomic_number() != new_element->atomic_number()))
  {
    _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;
    _implicit_hydrogens_known = 0;

    set_modified();
  }

  _permanent_aromatic = new_element->permanent_aromatic();

  _element = new_element;

  return;
}


angle_t
angle_between_atoms (const Atom & a1, const Atom & a2, const Atom & a3, const Atom & a4)
{
  assert (a1.ok());
  assert (a2.ok());
  assert (a3.ok());
  assert (a4.ok());

  Coordinates v21 (a1.x() - a2.x(), a1.y() - a2.y(), a1.z() - a2.z());
  Coordinates v32 (a2.x() - a3.x(), a2.y() - a3.y(), a2.z() - a3.z());
  Coordinates v43 (a3.x() - a4.x(), a3.y() - a4.y(), a3.z() - a4.z());

  v21.normalise();
  v32.normalise();
  v43.normalise();

  v21.cross_product (v32);
  v43.cross_product (v32);
  v43 *= static_cast<coord_t> (-1.0);

  angle_t rc = v21.angle_between (v43);

// Now we need to work out the directionality of the angle
// The cross product of v21 and v43 will be in the same or opposite
// direction of v32

  v21.cross_product(v43);

  angle_t tmp = v21.angle_between(v32);

  if (tmp < static_cast<angle_t>(0.1))
    return -rc;

  return rc;
}

/*
  Fully audit an atom. 
  Note the ok() function is inlined.
*/

int
Atom::audit() const
{
  if (0 != ok())
    return 1;

  if (! resizable_array<Bond *>::ok())
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (! _things[i]->ok())
      return 0;
  }

  return 1;
}

/*
  The molecule is being resized. we need to get rid of any
  bonds to atoms which will be disappearing
*/

int
Atom::molecule_being_resized (int new_size)
{
  int rc = 0;
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    const Bond * b = _things[i];
    if (b->a1() >= new_size || b->a2() >= new_size)
    {
      rc++;
      remove_item(i);
    }
  }

  if (rc)
    set_modified();

  return rc;
}

/*
  Given that THIS is atom number A, what is the other atom
  associated with the I'th bond
*/

atom_number_t
Atom::other (atom_number_t a, int i) const
{
  assert (ok());
  assert (ok_index(i));

  return _things[i]->other(a);
}

const Bond *
Atom::bond_to_atom (atom_number_t a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves(a))
      return _things[i];
  }

  return NULL;
}

bond_type_t
Atom::btype_to_connection (int i) const
{
  assert (ok());

  return _things[i]->btype();
}

bond_type_t
Atom::btype_to_atom (atom_number_t a) const
{
  assert (ok());

  const Bond * b = bond_to_atom(a);
  assert (b && b->ok());

  return b->btype();
}

int
Atom::other_and_type (atom_number_t my_atom_number,
                      int i, atom_number_t & a, bond_type_t & bt) const
{
  assert (ok());
  assert (ok_index (i));

  const Bond * b = _things[i];

  a  = b->other(my_atom_number);
  bt = b->btype();

  return 1;
}

int
Atom::connections (atom_number_t my_atom_number,
                   atom_number_t * others, bond_type_t * bt) const
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    others[i] = b->other(my_atom_number);

    if (NULL != bt)
      bt[i] = b->btype();
  }

  return _number_elements;
}

int
Atom::connections (atom_number_t my_atom_number,
                   Set_of_Atoms & others) const
{
  assert (ok());

  others.resize_keep_storage(0);
  others.make_room_for_extra_items(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];
    others.add(b->other(my_atom_number));
  }

  return _number_elements;
}

int
Atom::connections_and_types (atom_number_t my_atom_number,
                             Set_of_Atoms & others,
                             resizable_array<bond_type_t> & bt) const
{
  assert (ok());

  others.resize_keep_storage(0);
  others.make_room_for_extra_items(_number_elements);

  bt.resize_keep_storage(0);
  bt.make_room_for_extra_items(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];
    others.add(b->other(my_atom_number));
    bt.add(b->btype());
  }

  return _number_elements;
}

int
Atom::connections_and_types (atom_number_t my_atom_number,
                             atom_number_t * others,
                             bond_type_t * bt) const
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    others[i] = b->other(my_atom_number);
    bt[i] = b->btype();
  }

  return _number_elements;
}

int
Atom::bond_types (resizable_array<bond_type_t> & bt) const
{
  assert (ok());

  bt.resize_keep_storage(0);
  bt.resize_keep_storage(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    bt.add (_things[i]->btype());
  }

  return _number_elements;
}

int
Atom::bond_types (bond_type_t * bt) const
{
  assert (ok());
  for (int i = 0; i < _number_elements; i++)
  {
    bt[i] = _things[i]->btype();
  }

  return _number_elements;
}

Coordinates &
form_unit_vector (const Atom & a1, const Atom & a2)
{
  Coordinates * rc = new Coordinates(a1.x() - a2.x(),
                                     a1.y() - a2.y(),
                                     a1.z() - a2.z());
  rc->normalise();

  return * rc;
}

Coordinates &
form_vector (const Atom * a)
{
  assert (OK_ATOM (a));
  Coordinates * result = new Coordinates(a->x(), a->y(), a->z());

  return * result;
}

/*
  An atom is a resizable array of Connections
  We make the assumption that from any atom (THIS) there will be
  only one connection to any other atom (ZATOM)
*/

int
Atom::remove_bonds_to_atom (atom_number_t zatom, int adjust_atom_numbers)
{
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    Bond * b = _things[i];
    if (b->involves(zatom))
    {
      remove_item (i);
      if (0 == adjust_atom_numbers)
        break;
    }
    else if (adjust_atom_numbers)
      b->adjust_for_loss_of_atom(zatom);
  }

  _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;
  _implicit_hydrogens_known = 0;

  set_modified();

  return 0;     // atom not found bonded
}

int
Atom::is_isotope() const
{
  assert (ok());

  return _isotope;
}


int
Atom::is_bonded_to (atom_number_t a1) const
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves (a1))
      return 1;
  }

  return 0;
}

int
Atom::is_bonded_to (atom_number_t a1, bond_type_t & bt) const
{
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];
    if (b->involves(a1))
    {
      bt = b->btype();
      return 1;
    }
  }

  return 0;
}

/*
*/

void
Atom::set_modified()
{
  _nrings = ATOM_PROPERTY_UNKNOWN;

  if (! _implicit_hydrogens_known)
    _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;

  _nbonds = ATOM_PROPERTY_UNKNOWN;

  return;
}

int
Atom::set_formal_charge (formal_charge_t f)
{
  if (! reasonable_formal_charge_value(f))
  {
    if (_element->organic())
      cerr << "Atom::set_formal_charge: possibly unreasonable formal charge " << f << " on atom of type " << atomic_symbol() << endl;
  }

  _formal_charge = f;

  _implicit_hydrogens_known = 0;

  set_modified();

  return 1;
}

int
Atom::set_implicit_hydrogens (int h, int override_known_value)
{
  assert (h >= 0);

  if (_implicit_hydrogens_known && 0 == override_known_value)
  {
    cerr << "Atom::set_implicit_hydrogens: cannot override known value\n";
    iwabort();
    return 0;
  }

  _implicit_hydrogens = h;
  _implicit_hydrogens_known = 0;

  _nbonds = ATOM_PROPERTY_UNKNOWN;

  return 1;
}

void
Atom::set_implicit_hydrogens_known (int i)
{
  if (0 == i)
  {
    _implicit_hydrogens_known = 0;
    return;
  }

  if (ATOM_PROPERTY_UNKNOWN == _implicit_hydrogens)
    _implicit_hydrogens = 0;

  _implicit_hydrogens_known = i;

  _nbonds = ATOM_PROPERTY_UNKNOWN;

  return;
}

int
Atom::implicit_hydrogens()
{
  if (ATOM_PROPERTY_UNKNOWN == _implicit_hydrogens)
  {
    int tmp;
    (void) compute_implicit_hydrogens(tmp);
    _implicit_hydrogens = static_cast<short>(tmp);

    return tmp;     // save conversion from short to int below
  }

  return _implicit_hydrogens;
}

void
Atom::unset_all_implicit_hydrogen_information ()
{
  _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;
  _implicit_hydrogens_known = 0;

  return;
}

static int display_abnormal_valence_messages = 1;

void
set_display_abnormal_valence_messages (int d)
{
  display_abnormal_valence_messages = d;
}

static int
possible_hypervalent_hcount(const int v)
{
  if (alternate_valences_give_hcount)
  {
    assert (v >= 0);
    return v;
  }

  return 0;
}

//#define DEBUG_COMPUTE_IMPLICIT_HYDROGENS

int
Atom::compute_implicit_hydrogens (int & result)
{
  if (ATOM_PROPERTY_UNKNOWN == _nbonds)
    (void) nbonds();

  int ivalence = _element->normal_valence();

#ifdef DEBUG_COMPUTE_IMPLICIT_HYDROGENS
  cerr << "Computing implicit hydrogens for (" << _element->symbol();
  if (_formal_charge > 0)
    cerr << '+';
  else if (_formal_charge < 0)
    cerr << '-';
  cerr << ") ncon " << _number_elements << " bonds = " << _nbonds <<
          " valence = " << ivalence << endl;
#endif

  if (VALENCE_NOT_DEFINED == ivalence)
  {
    result = 0;
    return 1;
  }

  int adjusted_valence = ivalence + _formal_charge;

// Consider C+. It starts with 4 unpaired electrons. Removing an
// electron leaves just 3 unpaired electrons, whereas the computed
// adjusted valence above would be 5

  int ose;

  if (_element->outer_shell_electrons(ose) && adjusted_valence > ose)
  {
    int tmp = ose - _formal_charge - _nbonds;

    if (tmp > 3)     // The case of Cl=N-
    {
      if (display_abnormal_valence_messages)
        cerr << "Atom::compute_implicit_hydrogens: " << tmp << " Hydrogens, impossible " << _element->symbol() << " ncon " << _number_elements << endl;
      result = 0;
      return 1;
    }

    if (tmp >= 0)
    {
      result = static_cast<unsigned int>(tmp);
      return 1;
    }

    if (display_abnormal_valence_messages)
    {
      cerr << "Atom::compute_implicit_hydrogens: strange chemistry: " << _element->symbol() << endl;
      cerr << "Element has " << ose << " outer shell electrons, and normal valence of " <<
              _element->normal_valence() << endl;
      cerr << "Atom has " << _nbonds << " bonds";
      if (_formal_charge)
        cerr << ", and " << _formal_charge << " formal_charge";
      cerr << endl;
    }

    result = 0;
    return 0;
  }

  if (_formal_charge <= 0)
  {
    if (_nbonds - _formal_charge <= ivalence)
    {
      result = ivalence - (_nbonds - _formal_charge);
      return 1;
    }
  }
  else if (_nbonds - _formal_charge <= ivalence)
  {
    result = ivalence - (_nbonds - _formal_charge);
    return 1;
  }

#ifdef DEBUG_COMPUTE_IMPLICIT_HYDROGENS
  cerr << "Looking for alternate valence states\n";
#endif

// July 2002, ran into some Cl=N- molecules.

  if (17 == _element->atomic_number() && 1 == _number_elements && 0 == _formal_charge && 2 == _nbonds)
  {
    result = 0;    
    return 0;      // return 0 so it will be known as a valence error
  }

// 4 connected neutral Nitrogen is an optional behaviour

  if (7 == _element->atomic_number() && 4 == _number_elements && 0 == _formal_charge && 4 == nbonds())
  {
    result = file_scope_four_connected_neutral_nitrogen_has_h;
    return 1;
  }

// Look for an alternate valences

  int highest_available_valence = ivalence;

  for (int j = 0; j < _element->number_alternate_valences(); j++)
  {
    int ivalence = _element->alternate_valence(j);
    if (ivalence > highest_available_valence)
      highest_available_valence = ivalence;

    if (_formal_charge <= 0)
    {
      if (_nbonds <= ivalence - _formal_charge)
      {
        if (7 == _element->atomic_number() && 0 == _formal_charge && 4 == _nbonds && 3 == _number_elements)      // those horrible *-N(=*)-* molecules
          result = 0;
        else if (16 == _element->atomic_number() && -1 == _formal_charge && 2 == _number_elements)
          result = 1;
        else
          result = possible_hypervalent_hcount(ivalence - _formal_charge - _nbonds);
        return 1;
      }
    }
    else if (ivalence - _formal_charge - _nbonds >= 0)
    {
      result = possible_hypervalent_hcount(ivalence - _formal_charge - _nbonds);
      return 1;
    }
  }

// This next case comes from seeing a [P-]F6. Jan 2010, now handled properly above

/*if (_formal_charge < 0 && highest_available_valence - _formal_charge == _nbonds)
  {
    result = 0;
    return 1;
  }*/

  if (! _element->organic())
  {
    result = 0;
    return 1;
  }

// Nov 98. the case -C#O fails all the other tests.

  if (8 == _element->atomic_number() && 1 == _number_elements && 3 == _nbonds)
  {
    result = 0;
    return 1;
  }

// If we come to here, the atom truly does seem to have an abnormally high valence.

  if (display_abnormal_valence_messages)
  {
    cerr << "Atom '" << _element->symbol() << "', has " << _nbonds << " bonds";
    if (_formal_charge)
      cerr << ", and charge " << _formal_charge;
    if (_element->number_alternate_valences())
      cerr << ", but highest valence is " << highest_available_valence << endl;
    else
      cerr << ", but valence " << ivalence << endl;
  }

  result = 0;
  return 0;
}

void
Atom::invalidate_computed_values_after_bond_change()
{
  _nbonds = ATOM_PROPERTY_UNKNOWN;

  if (! _implicit_hydrogens_known)
    _implicit_hydrogens = ATOM_PROPERTY_UNKNOWN;

  return;
}

int
Atom::recompute_implicit_hydrogens (int & result)
{
  if (! compute_implicit_hydrogens(result))
    return 0;

  _implicit_hydrogens = result;

  return 1;
}

int
Atom::lone_pair_count (int & result)
{
  int electrons;
  if (! _element->outer_shell_electrons(electrons))
    return 0;

  electrons -= _formal_charge;
  assert (electrons >= 0);

  int rc = electrons - nbonds() - implicit_hydrogens();
  if (rc < 0)
  {
    cerr << "Atom::lone_pair_count: Yipes, electron count error, type " << _element->symbol() << endl;

    (void) _element->outer_shell_electrons(electrons);
    cerr << "Outer shell electrons is " << electrons << endl;

    cerr << "Formal charge is " << _formal_charge << endl;
    cerr << "Nbonds is " << nbonds() << endl;
    cerr << "Implicit hydrogens is " << _implicit_hydrogens << endl;
    cerr << "Rc = " << rc << endl;
    cerr << "Likely invalid valence\n";
    return 0;
  }

// Interesting point here. If it has one unpaired electron, should that
// count as a lone pair. Similarly if it has three unpaired electrons,
// is that one or two lone pairs. Let's try it this way and see what happens

  result = rc / 2;

// 3 connected iodine must be a special case (NCI dat aug 97, got
// several chiral, 3 connected Iodine atoms

// Jul 02. My old friend 3 connected neutral Nitrogen with 4 bonds. Give it a lone pair

// Sep 02. Dan Robertson had some 4 connected neutral Nitrogen atoms

  atomic_number_t z = _element->atomic_number();

  if (6 == z)     // intercept the most common case here for efficiency
    ;
  else if (53 == z && 3 == _number_elements && 3 == nbonds() && 2 == result)
    result = 1;
  else if (7 == z && 3 == _number_elements && 0 == _formal_charge && 4 == nbonds())
    result = 1;
  else if (7 == z && 4 == _number_elements && 0 == _formal_charge && 4 == nbonds())
    result = 0;

  return 1;
}

int
Atom::pi_electrons (int & pe)
{
  if (ATOM_PROPERTY_UNKNOWN == _implicit_hydrogens)
    (void) implicit_hydrogens();

  return _element->pi_electrons(_number_elements + _implicit_hydrogens,   // total number of connections
                                _formal_charge, pe);
}

//#define DEBUG_ATOM_VALENCE_OK

int
Atom::valence_ok()
{
  if (ATOM_PROPERTY_UNKNOWN == _implicit_hydrogens)
    (void) implicit_hydrogens();

#ifdef DEBUG_ATOM_VALENCE_OK
  cerr << "Checking valence of atom " << _element->symbol() << endl;
#endif

  int ih;
  if (! compute_implicit_hydrogens(ih))
    return 0;

#ifdef DEBUG_ATOM_VALENCE_OK
    cerr << "IH " << ih << " ncon " << _number_elements << endl;
#endif
  if (_number_elements + ih > 4)
  {
    if (15 == _element->atomic_number() && 6 == (_number_elements + ih))
      ;
    else
      return 0;
  }

  if (! reasonable_formal_charge_value(_formal_charge))
    return 0;

  if (! _implicit_hydrogens_known)
    _implicit_hydrogens = ih;

// Mar 99. Need to be more careful than previously about this. We need to
// check whether or not the connections present satisfy one of the known
// valence states.

  int nb = nbonds() + _implicit_hydrogens - _formal_charge;

#ifdef DEBUG_ATOM_VALENCE_OK
  cerr << "z = " << _element->atomic_number() << " final check, " << _implicit_hydrogens << " implicit hydrogens, " << nbonds() << " bonds. fc = " << _formal_charge << ". nb = " << nb << endl;
#endif

  if (_element->is_valid_valence(nb))
    return 1;

// We don't check the valence of strange elements

  if (! _element->organic())
    return 1;

  int ose;
  if (! _element->outer_shell_electrons(ose))
    return 0;

#ifdef DEBUG_ATOM_VALENCE_OK
  cerr << "Atom has " << ose << " outer shell electrons\n";
#endif

// consider the case of 3 connected C+. nb = 2

  if (_formal_charge > 0)
  {
    if (nbonds() + _implicit_hydrogens == ose - _formal_charge)    // atom now has lone pair
      return 1;
  }

// Perchlorates will be entered as Cl- with 4 =O around them

  if (8 == nbonds() && (8 == ose - _formal_charge))
    return 1;

// Drive out positively charged halogens

  if (_element->is_halogen() && _formal_charge > 0)
    return 0;

// Get PF6-

  if (6 == _number_elements && 6 == nbonds() && -1 == _formal_charge && 0 == ih)
    return 1;

  return 0;
}

int
Atom::in_another_ring()
{
  if (ATOM_PROPERTY_UNKNOWN == _nrings)
    _nrings = 1;
  else
    _nrings++;

  return _nrings;
}

static int
bond_order_comparitor (Bond * const * pb1, Bond * const * pb2)
{
  const Bond * b1 = *pb1;
  const Bond * b2 = *pb2;

  if (b1->btype() == b2->btype())
    return 0;

  if (b1->is_aromatic())
  {
    if (b2->is_aromatic())
      return 0;

    return -1;
  }
  else if (b2->is_aromatic())
    return 1;


  int nb1 = b1->number_of_bonds();
  int nb2 = b2->number_of_bonds();

  if (nb1 < nb2)
    return -1;
  else if (nb1 == nb2)
    return 0;
  else
    return 1;
}

void
Atom::multiple_bonds_first()
{
  if (_number_elements <= 1)
    return;

  if (2 == _number_elements)
  {
    Bond * b0 = _things[0];
    const Bond * b1 = _things[1];

    if (b0->btype() == b1->btype() || (b0->is_aromatic() && b1->is_aromatic()))
      return;

    if (b1->is_single_bond() || b1->is_aromatic())
      return;

    int nb0 = b0->number_of_bonds();
    int nb1 = b1->number_of_bonds();

    if (nb0 < nb1)
    {
      _things[0] = _things[1];
      _things[1] = b0;
    }

    return;
  }

  sort (bond_order_comparitor);

  return;
}

int
Atom::next_atom_for_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const
{
  (void) canonical_order;    // not used

  atom_number_t zdefault = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    atom_number_t j = b->other(my_atom_number);
    if (already_done[j])
      continue;

    if (! b->is_single_bond())
    {
      next_atom = j;
      return 1;
    }
    else if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;
  }

  if (INVALID_ATOM_NUMBER == zdefault)
    return 0;

  next_atom = zdefault;
  return 1;
}

int
Atom::next_atom_for_unique_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const
{
  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int rsave = 0;                 // initialised to shut gcc up

  int highest_bond_count = 0;     // initialised to shut gcc up

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    atom_number_t j = b->other(my_atom_number);

    if (already_done[j] >= 0)
      continue;

    int bcount;

//  Note that the counting here is different from in the function number_of_bonds().
//  In that function, single and double bonds are counted first, here we want to catch
//  aromatic bonds first as contributing one to the bond count

    if (b->is_aromatic())
      bcount = 1;
    else if (b->is_single_bond())
      bcount = 1;
    else if (b->is_double_bond())
      bcount = 2;
    else if (b->is_triple_bond())
      bcount = 3;
    else         // GACK!!
      bcount = (int) b->btype();

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Choosing next unique atom, " << b << " bond " << bt << " (bcount = " << bcount << ", rb = " << rb << ")\n";
#endif

    int rj = canonical_order[j];

    if (INVALID_ATOM_NUMBER == zdefault || bcount > highest_bond_count)
    {
      highest_bond_count = bcount;
      zdefault = j;
      rsave = rj;
    }
    else if (highest_bond_count == bcount && rj > rsave)
    {
#ifdef DEBUG_UNIQUE_SMILES_ORDERING
      cerr << "Bonds equal, rj = " << rn << " rsave = " << rsave << endl;
#endif

      zdefault = j;
      rsave = rj;
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Next unique atom is atom " << next_atom << endl;
#endif

    return 1;
  }

  return 0;
}

static Random_Number_Working_Storage smiles_random_number_stream;

int
Atom::next_atom_for_random_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const
{
  (void) canonical_order;    // not used

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  atom_number_t multiple_bond = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    atom_number_t j = b->other(my_atom_number);
    if (already_done[j] >= 0)
      continue;

    if (! b->is_single_bond())
    {
      if (INVALID_ATOM_NUMBER == multiple_bond)
        multiple_bond = j;
      else if (smiles_random_number_stream.random_one_or_zero())
        multiple_bond = j;
    }
    else if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;
    else if (smiles_random_number_stream.random_one_or_zero())
      zdefault = j;
  }

  if (INVALID_ATOM_NUMBER != multiple_bond)
  {
    next_atom = multiple_bond;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;
    return 1;
  }

  return 0;
}

int
Atom::set_bond_type_to_atom (atom_number_t zatom, bond_type_t bt)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves(zatom))
    {

      _things[i]->set_bond_type(bt);

      set_modified();

      return 1;
    }
  }

  cerr << "Atom::set_bond_type_to_atom: atom not bonded to atom " << zatom << endl;

  assert (NULL == "This is very bad");

  return 0;
}

int
Atom::is_halogen() const
{
  return _element->is_halogen();
}

/*
  Could not get exactly formatted I/O with the SGI C++ library
  Lots of problems with sprintf in the Sun run-time library too
*/

int
Atom::write_coordinates (ostream & os, 
                         int include_space) const
{
#if defined(__GNUG__) || defined (IW_INTEL_COMPILER)
  ios::fmtflags old_flags = os.flags(ios::fixed);
#else
  long old_flags = os.flags(ios::fixed);
#endif

  int old_precision  = os.precision(4);

  if (include_space)
    os << setw(10) << _x << ' ' << setw(10) << _y << ' ';
  else
    os << setw(10) << _x << setw(10) << _y;

// Writing 0.0 Z coordinates is common when writing 2d files

  if (static_cast<coord_t>(0.0) == _z)
    os << "    0.0000 ";
  else
    os << setw(10) << _z << ' ';
  
  os.precision(old_precision);
  os.flags(old_flags);

/*#ifdef sun
  char buffer[80];

  IW_SPRINTF (buffer, "%10.4f", _x);
  os << buffer;
  IW_SPRINTF (buffer, "%10.4f", _y);
  os << buffer;
  IW_SPRINTF (buffer, "%10.4f", _z);
  os << buffer << ' ';

#else
  IW_SPRINTF (buffer, "%10.4f%10.4f%10.4f ", _x, _y, _z);
  os << buffer;
#endif*/

  return os.good();
}

int
Atom::ncon (atom_number_t zatom, const int * include_atom) const
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];
    
    atom_number_t j = b->other(zatom);

    if (include_atom[j])
      rc++;
  }

  return rc;
}

int
Atom::is_centre_of_cis_trans_bond() const
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Bond * b = _things[i];

    if (b->is_double_bond() && b->part_of_cis_trans_grouping())
      return 1;
  }

  return 0;
}

int
Atom::number_directional_bonds_attached() const
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->is_directional())
      rc++;
  }

  return rc;
}
void reset_atom_file_scope_variables()
{
  current_atoms=0;
  copy_implicit_hydrogen_count_in_atom_copy_constructor = 1;
  formal_charge_t _min_reasonble_atomic_formal_charge_value = -5;
  formal_charge_t _max_reasonble_atomic_formal_charge_value =  5;
}
