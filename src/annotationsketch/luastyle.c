/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "lua.h"
#include "annotationsketch/luastyle.h"

/* key used to store the Style object in the Lua registry */
#define STYLE_KEY gt_style_new

void lua_put_style_in_registry(lua_State *L, GtStyle *style)
{
  assert(L && style);
  lua_pushlightuserdata(L, STYLE_KEY);
  lua_pushlightuserdata(L, style);
  lua_rawset(L, LUA_REGISTRYINDEX);
}

GtStyle* lua_get_style_from_registry(lua_State *L)
{
  GtStyle *style;
  assert(L);
  lua_pushlightuserdata(L, STYLE_KEY);
  lua_rawget(L, LUA_REGISTRYINDEX);
  assert(lua_islightuserdata(L, -1));
  style = lua_touserdata(L, -1);
  lua_pop(L, 1);
  return style;
}
