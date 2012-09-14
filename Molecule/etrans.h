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
#ifndef IW_ELEMENT_TRANS_H
#define IW_ELEMENT_TRANS_H

#include <iostream>

#include "iwaray.h"

class Element;
class Molecule;
class IWString;
class Molecule_to_Match;

#include "ematch.h"

class Element_Transformation
{
  private:
    int   _transform_every_atom_type;
    Element_Matcher _from;
    const Element * _to;
    int _isotope;

    int _molecules_processed;
    int _molecules_changed;
    int _atoms_changed;

//  private functions

    void _default_values ();

  public:
    Element_Transformation ();

    int ok () const;
    int debug_print (ostream &) const;

//  int build (const char *);
    int build (const IWString &);

    int process (Molecule &);

    int process (Molecule_to_Match &);
};

class Element_Transformations : public resizable_array_p<Element_Transformation>
{
  private:
  public:

    int ok () const;
    int debug_print (ostream &) const;

    int active () const { return _number_elements;}

    int construct_from_command_line (Command_Line &, int = 0, char = 't');

    int process (Molecule *);

    int process (Molecule &);

    int process (Molecule_to_Match &);
};

class Command_Line;

extern int display_standard_etrans_options (ostream &, char = 't');

extern int process_element_transformations (Command_Line &,
                                            Element_Transformations &,
                                            int = 0,
                                            char = 't');
#endif
