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
  Strips leading spaces
*/

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"

void
strip_leading_blanks (char *cc)
{
  char *c = cc;
  while (isspace (*c) && '\0' != *c)
    c++;

  if (c == cc)    /* no whitespace found */
    return;
  while (*cc)
  {
    *cc++ = *c++;
  }
  *cc = '\0';
}

void
strip_trailing_blanks (char *cc)
{
  size_t i = strlen (cc);
  char *c = cc + i - 1;

  if (0 == i)
    return;

  while (isspace (*c))
  {
    if (c == cc)
    {
      *c = '\0';
      return;
    }
    c--;
  }
  
  c++;
  *c = '\0';
  return;
}

void
no_newline (char *cc)
{
  char *c = strchr (cc, '\n');
  if (NULL != c)
    *c = '\0';
  return;
}

char *
basename (char *path_name)
{
  char *c, *file_name;

  c = file_name = path_name;
  while (*c)
  {
    if (*c == '/')
      file_name = c + 1;
    c++;
  }

  return (file_name);
}

int
remove_blanks (char *cc)
{
  char *c1, *c2;
  int nblanks = 0;

  c1 = c2 = cc;

  while (*c1)
  {
    if (' ' == *c1)
    {
      nblanks++;
      c1++;
    }
    else
    {
      *c2 = *c1;
      c1++;
      c2++;
    }
  }
  *c2 = '\0';
  return nblanks;
}

void
to_lowercase (char *cc, int characters_to_convert)
{
  int characters_processed = 0;
  while (*cc)
  {
    *cc = tolower (*cc);
    cc++;
    characters_processed++;
    if (characters_to_convert && characters_processed >= characters_to_convert)
      return;
  }

  return;
}

/*void
to_lowercase (resizable_array_p<char> *buffers)
{
  for (int i = 0; i < buffers->number_elements (); i++)
    to_lowercase (buffers->item (i));
}*/

void
to_uppercase (char *cc, int characters_to_convert)
{
  int characters_processed = 0;
  while (*cc)
  {
    *cc = toupper (*cc);
    cc++;
    characters_processed++;
    if (characters_to_convert > 0 && characters_processed >= characters_to_convert)
      return;
  }

  return;
}

/*void
to_uppercase (resizable_array_p<char> *buffers)
{
  for (int i = 0; i < buffers->number_elements (); i++)
    to_uppercase (buffers->item (i));
}*/

int
compress_blanks (char * zstring)
{
  int from = 0;
  int to = 0;

  int rc = 0;
  int previous_was_space = 0;
  while (zstring[from])
  {
    int copy_char = 1;

    if (isspace (zstring[from]))
    {
      if (previous_was_space)
      {
        copy_char = 0;
        rc++;
      }
      previous_was_space = 1;
    }
    else
      previous_was_space = 0;

    if (copy_char)
    {
      zstring[to] = zstring[from];
      to++;
    }

    from++;
  }

  zstring[to] = '\0';

  return rc;
}

/*
  We return the length of the string
*/

int
compress_blanks (char * zstring, int nchars)
{
  int from = 0;
  int to = 0;

  int previous_was_space = 0;
  for (int i = 0; i < nchars; i++)
  {
    int copy_char = 1;

    if (isspace (zstring[from]))
    {
      if (previous_was_space)
        copy_char = 0;
      previous_was_space = 1;
    }
    else
      previous_was_space = 0;

    if (copy_char)
    {
      zstring[to] = zstring[from];
      to++;
    }

    from++;
  }

  return to;
}

int
is_int (const char *buffer, int *i)
{
  char *c;
  int tmp = strtol (buffer, &c, 10);

  if (c == buffer)
    return 0;

  if ('\0' == *c)
    ;     // good
  else if (! isspace (*c))
    return 0;

  *i = tmp;
  return 1;
}

int
is_double (const char *buffer, double *x)
{
  assert (NULL != buffer);
  assert (NULL != x);
  char *c;
  double tmp = strtod (buffer, &c);

//cerr << "Parsing '" << buffer << "' as double\n";

  if (c == buffer)
    return 0;

// c must be pointing to a string terminator, or whitespace

  if ('\0' == *c)
    ;     // good
  else if (! isspace (*c))
    return 0;

  *x = tmp;
  return 1;
}

/*
  Counts the number of occurrences of NEEDLE in HAYSTACK
*/

int
ccount (const char *haystack, char needle)
{
  assert (NULL != haystack);

  int count = 0;

  while (*haystack) 
  {
    if (needle == *haystack)
      count++; 
    haystack++;
  }
    
  return count;
} 

char *
make_copy (const char *c)
{
  assert (NULL != c);
  int lenc = static_cast<int>(strlen (c));

  char * copyc = new char[lenc + 1];
  memcpy (copyc, c, lenc);

  return copyc;
}

int
words (const char *string, const char separator)
{

  int in_word;

  if ('\0' == *string)
    return 0;        // null string has no words

  if (*string == separator)
    in_word = 0;
  else
    in_word = 1;

  int words = 0;

  while (1)
  {
    string++;

    if ('\0' == *string)
    {
      if (in_word)
        words++;
      return words;
    }

    if (separator == *string)
    {
      if (in_word)
      {
        words++;
        in_word = 0;
      }
    }
    else
      in_word = 1;
  }
}

const char *
iwbasename (const char * fname)
{
  const char * rc = fname;
  int previous_char_was_slash = 0;

  while (*fname)
  {
    if ('/' == *fname)
      previous_char_was_slash = 1;
    else
    {
      if (previous_char_was_slash)
        rc = fname;
      previous_char_was_slash = 0;
    }
    fname++;
  }

  return rc;
}
