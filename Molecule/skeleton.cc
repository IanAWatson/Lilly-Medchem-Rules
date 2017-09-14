/*
  Compute various things about conjugated sections in a molecule.
*/

#include <stdlib.h>
#include <memory>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Message about what the programme does\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
your_programme (Molecule & m,
                IWString_and_File_Descriptor & output)
{
  int nitrogen_count = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 == m.atomic_number(i))
      nitrogen_count++;
  }

  output << m.name() << " has " << nitrogen_count << " nitrogen atoms\n";

  return output.good();
}

static int
your_programme (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! your_programme(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
your_programme (const char * fname, int input_type, 
                IWString_and_File_Descriptor & output)
{
  assert(NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return your_programme(input, output);
}

static int
your_programme (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:l");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! your_programme(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = your_programme(argc, argv);

  return rc;
}
