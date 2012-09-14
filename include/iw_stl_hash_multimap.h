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

#if (__GNUC__ >= 3)
#include <ext/hash_map>
using namespace __gnu_cxx;
#else
#include <hash_map>
using namespace stdext;
#endif

using namespace std;

#include "iwstring.h"
#include "iwhash.h"

template <typename K, typename V>
class IW_STL_Hash_Multimap : public hash_multimap<K, V, IWStringHash>
{
  private:
  public:
};

typedef IW_STL_Hash_Multimap<IWString, int> IW_STL_Hash_Multimap_int;
typedef IW_STL_Hash_Multimap<IWString, float> IW_STL_Hash_Multimap_float;
typedef IW_STL_Hash_Multimap<IWString, double> IW_STL_Hash_Multimap_double;
typedef IW_STL_Hash_Multimap<IWString, IWString> IW_STL_Hash_Multimap_IWString;

#endif
