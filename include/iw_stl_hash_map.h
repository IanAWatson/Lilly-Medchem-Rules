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
#ifndef IW_STL_HASH_MAP_H
#define  IW_STL_HASH_MAP_H

#include <stdlib.h>

#include "iwconfig.h"

#if (GCC_VERSION >= 40400)
#include <unordered_map>
#define IW_Hash_Map std::unordered_map
#elif (__GNUC__ >= 3)
#include <ext/hash_map>
using namespace __gnu_cxx;
#define IW_Hash_Map hash_map
#else
#include <hash_map>
using namespace stdext;
#define IW_Hash_Map hash_map
#endif

using namespace std;

#include "iwstring.h"
#include "iwhash.h"

template <typename T, typename V>
class IW_STL_Hash_Map : public IW_Hash_Map<T, V, IWStringHash>
{
  private:
  public:
    int contains (const T & t) const 
    {
      typename IW_Hash_Map<T, V, IWStringHash>::const_iterator f = IW_Hash_Map<T, V, IWStringHash>::find (t);

      return f != IW_Hash_Map<T, V, IWStringHash>::end ();
    }

    typedef typename IW_STL_Hash_Map<T, V>::const_iterator const_iterator;
};

typedef IW_STL_Hash_Map<IWString, int> IW_STL_Hash_Map_int;
typedef IW_STL_Hash_Map<IWString, unsigned int> IW_STL_Hash_Map_uint;
typedef IW_STL_Hash_Map<IWString, float> IW_STL_Hash_Map_float;
typedef IW_STL_Hash_Map<IWString, long> IW_STL_Hash_Map_long;
typedef IW_STL_Hash_Map<IWString, off_t> IW_STL_Hash_Map_off_t;
typedef IW_STL_Hash_Map<IWString, IWString> IW_STL_Hash_Map_String;

#endif
