#ifndef IW_RX_LREGEX_H
#define IW_RX_LREGEX_H

#include <regex.h>

#include "iwstring.h"

/*
  This class is used with IW_Regular_Expression_Template
*/

class IW_lregex
{
  private:
    regex_t  _compiled_regular_expression;
    int      _regcomp_flags;
    int      _reg_eflags;

//  We need to keep track of whether or not _compiled_regular_expression
//  contains a valid compiled regular expression

    int _compiled;

 // We need to keep track of whether or not the previous match () filled
 // the _regmatch array

    int _regmatch_valid;

//  When asked to save the subexpressions, we also need to save
//  the string which was matched. Note that we store it as a const_IWSubstring
//  which means that the string must be kept in scope by the caller.

    regmatch_t * _regmatch;
    const_IWSubstring _string_which_was_matched;

//  private functions

    void _default_values ();

    int _matches_save_subexpressions (const char *, int);

  public:
    IW_lregex ();
    IW_lregex (const char *);
    IW_lregex (const char *, int);
    IW_lregex (const const_IWSubstring &);
    IW_lregex (const IWString &);
    ~IW_lregex ();

    int ok () const;

    int set_pattern (const char *, int);
    int set_pattern (const char *);
    int set_pattern (const IWString &);
    int set_pattern (const const_IWSubstring &);

    int matches (const char *);
    int matches (const char *, int);
    int matches (const IWString &);
    int matches (const const_IWSubstring &);

    int matches_save_subexpressions (const char *);
    int matches_save_subexpressions (const char *, int);
    int matches_save_subexpressions (const IWString &);
    int matches_save_subexpressions (const const_IWSubstring &);

    const regmatch_t * regmatch () const;
    int number_subexpressions () const;

    int dollar (int whichdollar, const_IWSubstring & zresult) const;

    const const_IWSubstring dollar (int) const;

    int dollar (int whichdollar, int &);
    int dollar (int whichdollar, float &);
    int dollar (int whichdollar, double &);

    int s (IWString &, const const_IWSubstring &);
    int s (const IWString &, const const_IWSubstring &, IWString &);

    void set_use_extended (int);
};

// In place conversion

extern int iwcrex_s (IWString &, IW_lregex &, const const_IWSubstring &);

// Result put in a new variable

extern int iwcrex_s (const IWString &, IW_lregex &, const const_IWSubstring &, IWString &);

// The Regexp is built on the fly. Hopefully this won't be called multiple times,
// as it must be compiled each time.

extern int iwcrex_s (IWString &, const const_IWSubstring &, const const_IWSubstring &);

#endif
