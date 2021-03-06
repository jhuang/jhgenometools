/*
  Copyright (c) 2007-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "annotationsketch/color_api.h"
#include "annotationsketch/style.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/thread.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/luahelper.h"
#include "extended/luaserialize.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/gt_lua.h"

struct GtStyle
{
  lua_State *L;
  unsigned long reference_count;
  GtRWLock *lock, *clone_lock;
  bool unsafe;
  char *filename;
};

static inline void style_lua_new_table(lua_State *L, const char *key)
{
  lua_pushstring(L, key);
  lua_newtable(L);
  lua_settable(L, -3);
}

static const luaL_Reg luasecurelibs[] = {
  /* Think very hard before adding additional Lua libraries to this list, it
     might compromise the security of web applications like GenomeViewer!
     Do not add the 'io', 'os', or 'debug' library under any circumstances!
     Use the luainsecurelibs list for that. */
  {"", luaopen_base},
  {LUA_TABLIBNAME, luaopen_table},
  {LUA_STRLIBNAME, luaopen_string},
  {LUA_MATHLIBNAME, luaopen_math},
  {"gt", gt_lua_open_lib}, /* open the GenomeTools library for callbacks */
  {NULL, NULL}
};

static const luaL_Reg luainsecurelibs[] = {
  /* These are functions affecting the system outside the Lua sandbox!
     They should only be loaded in environments where unwanted code execution
     is not a security issue! */
  {LUA_OSLIBNAME, luaopen_os},
  {LUA_IOLIBNAME, luaopen_io},
  {NULL, NULL}
};

static void luaL_opencustomlibs(lua_State *L, const luaL_Reg *lib)
{
  for (; lib->func; lib++) {
    lua_pushcfunction(L, lib->func);
    lua_pushstring(L, lib->name);
    lua_call(L, 1, 0);
  }
}

GtStyle* gt_style_new(GtError *err)
{
  GtStyle *sty;
  gt_error_check(err);
  sty = gt_calloc(1, sizeof (GtStyle));
  sty->filename = NULL;
  sty->L = luaL_newstate();
  sty->lock = gt_rwlock_new();
  sty->unsafe = false;
  sty->clone_lock = gt_rwlock_new();
  if (!sty->L) {
    gt_error_set(err, "out of memory (cannot create new Lua state)");
    gt_free(sty);
    return NULL;
  }
  else
    luaL_opencustomlibs(sty->L, luasecurelibs);
  return sty;
}

GtStyle* gt_style_new_with_state(lua_State *L)
{
  GtStyle *sty;
  gt_assert(L && !lua_gettop(L)); /* make sure the Lua stack is empty */
  sty = gt_calloc(1, sizeof (GtStyle));
  sty->L = L;
  sty->unsafe = true;
  sty->lock = gt_rwlock_new();
  return sty;
}

GtStyle* gt_style_ref(GtStyle *style)
{
  gt_assert(style);
  gt_rwlock_wrlock(style->lock);
  style->reference_count++;
  gt_rwlock_unlock(style->lock);
  return style;
}

void gt_style_unsafe_mode(GtStyle *style)
{
  gt_assert(style);
  gt_rwlock_wrlock(style->lock);
  luaL_opencustomlibs(style->L, luainsecurelibs);
  style->unsafe = true;
  gt_rwlock_unlock(style->lock);
}

void gt_style_safe_mode(GtStyle *style)
{
#ifndef NDEBUG
  int stack_size;
#endif
  const luaL_Reg *lib = luainsecurelibs;
  gt_assert(style);
  gt_rwlock_wrlock(style->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(style->L);
#endif
  for (; lib->name; lib++) {
    lua_pushnil(style->L);
    lua_setglobal(style->L, lib->name);
  }
  style->unsafe = false;
  gt_assert(lua_gettop(style->L) == stack_size);
  gt_rwlock_unlock(style->lock);
}

bool gt_style_is_unsafe(GtStyle *sty)
{
#ifndef NDEBUG
  int stack_size;
#endif
  const luaL_Reg *lib = luainsecurelibs;
  bool safe = true;
  gt_assert(sty);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  for (; safe && lib->name; lib++) {
    lua_getglobal(sty->L, lib->name);
    if (!lua_isnil(sty->L, -1)) {
      safe = false;
    }
    lua_pop(sty->L, 1);
  }
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return !safe;
}

int gt_style_load_file(GtStyle *sty, const char *filename, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int had_err = 0;
  gt_error_check(err);
  gt_rwlock_wrlock(sty->lock);
  gt_assert(sty && sty->L && filename);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  gt_rwlock_unlock(sty->lock);
  gt_rwlock_wrlock(sty->lock);
  sty->filename = gt_cstr_dup(filename);
  gt_log_log("Trying to load style file: %s...", filename);
  if (luaL_loadfile(sty->L, filename) || lua_pcall(sty->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run style file: %s", lua_tostring(sty->L, -1));
    had_err = -1;
    lua_pop(sty->L, 1);
  }
  if (!had_err) {
    lua_getglobal(sty->L, "style");
    if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
      gt_error_set(err, "'style' is not defined or is not a table in \"%s\"",
                   filename);
      had_err = -1;
    }
    lua_pop(sty->L, 1);
  }
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return had_err;
}

void gt_style_reload(GtStyle *sty)
{
  int rval;
  gt_assert(sty && sty->filename);
  rval = gt_style_load_file(sty, sty->filename, NULL);
  gt_assert(!rval); /* should not happen, file was loaded before */
}

/* Searches for <section> inside the style table, creating it if it does not
   exist and finally pushing it on the Lua stack (at the top).
   Returns the total number of items pushed on the stack by this function. */
static inline int style_find_section_for_setting(GtStyle* sty,
                                                 const char *section)
{
  int depth = 0;
  gt_assert(sty && section);
  lua_getglobal(sty->L, "style");
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    lua_newtable(sty->L);
    lua_setglobal(sty->L, "style");
    lua_getglobal(sty->L, "style");
  }
  depth++;
  lua_getfield(sty->L, -1, section);
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    style_lua_new_table(sty->L, section);
    lua_getfield(sty->L, -1, section);
  }
  depth++;
  return depth;
}

/* Searches for <section> inside the style table, returning -1 if it is not
   found. Otherwise the number of items pushed onto the stack is returned. */
static inline int style_find_section_for_getting(const GtStyle *sty,
                                                 const char *section)
{
  int depth = 0;
  gt_assert(sty && section);
  lua_getglobal(sty->L, "style");
  if (lua_isnil(sty->L, -1)) {
    gt_log_log("'style' is not defined");
    lua_pop(sty->L, 1);
    return -1;
  } else depth++;
  lua_getfield(sty->L, -1, section);
  if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
    lua_pop(sty->L, 2);
    return -1;
  } else depth++;
  return depth;
}

GtStyleQueryStatus gt_style_get_color(const GtStyle *sty, const char *section,
                                      const char *key, GtColor *color,
                                      GtFeatureNode *gn, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && section && key && color);
  gt_error_check(err);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  /* set default colors */
  color->red = 0.5; color->green = 0.5; color->blue = 0.5; color->alpha = 0.5;
  /* get section */
  i = style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  }
  /* lookup color entry for given feature */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1))
  {
    if (gn) {
      GtGenomeNode *gn_lua = gt_genome_node_ref((GtGenomeNode*) gn);
      gt_lua_genome_node_push(sty->L, gn_lua);
    }
    if (lua_pcall(sty->L, gn ? 1 : 0, 1, 0) != 0)
    {
      gt_error_set(err, "%s", lua_tostring(sty->L, -1));
      lua_pop(sty->L, 3);
      gt_assert(lua_gettop(sty->L) == stack_size);
      gt_rwlock_unlock(sty->lock);
      return GT_STYLE_QUERY_ERROR;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
    lua_pop(sty->L, 3);
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  } else i++;
  /* update color struct */
  lua_getfield(sty->L, -1, "red");
  if (!lua_isnil(sty->L, -1) && lua_isnumber(sty->L, -1))
    color->red = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  lua_getfield(sty->L, -1, "green");
  if (!lua_isnil(sty->L, -1) && lua_isnumber(sty->L, -1))
    color->green = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  lua_getfield(sty->L, -1, "blue");
  if (!lua_isnil(sty->L, -1) && lua_isnumber(sty->L, -1))
    color->blue = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  lua_getfield(sty->L, -1, "alpha");
  if (!lua_isnil(sty->L, -1) && lua_isnumber(sty->L, -1))
    color->alpha = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return GT_STYLE_QUERY_OK;
}

void gt_style_set_color(GtStyle *sty, const char *section, const char *key,
                        const GtColor *color)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && section && key && color);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  i = style_find_section_for_setting(sty, section);
  lua_getfield(sty->L, -1, key);
  i++;
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    style_lua_new_table(sty->L, key);
    lua_getfield(sty->L, -1, key);
  }
  lua_pushstring(sty->L, "red");
  lua_pushnumber(sty->L, color->red);
  lua_settable(sty->L, -3);
  lua_pushstring(sty->L, "green");
  lua_pushnumber(sty->L, color->green);
  lua_settable(sty->L, -3);
  lua_pushstring(sty->L, "blue");
  lua_pushnumber(sty->L, color->blue);
  lua_settable(sty->L, -3);
  lua_pushstring(sty->L, "alpha");
  lua_pushnumber(sty->L, color->alpha);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
}

GtStyleQueryStatus gt_style_get_str(const GtStyle *sty, const char *section,
                                    const char *key, GtStr *text,
                                    GtFeatureNode *gn, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && key && section);
  gt_error_check(err);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  /* get section */
  i = style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1))
  {
    if (gn) {
      GtGenomeNode *gn_lua = gt_genome_node_ref((GtGenomeNode*) gn);
      gt_lua_genome_node_push(sty->L, gn_lua);
    }
    if (lua_pcall(sty->L, gn ? 1 : 0, 1, 0) != 0)
    {
      gt_error_set(err, "%s", lua_tostring(sty->L, -1));
      lua_pop(sty->L, 3);
      gt_assert(lua_gettop(sty->L) == stack_size);
      gt_rwlock_unlock(sty->lock);
      return GT_STYLE_QUERY_ERROR;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_isstring(sty->L, -1)) {
    lua_pop(sty->L, i+1);
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  } else i++;
  /* retrieve string */
  gt_str_set(text, lua_tostring(sty->L, -1));
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return GT_STYLE_QUERY_OK;
}

void gt_style_set_str(GtStyle *sty, const char *section, const char *key,
                      GtStr *value)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && section && key && value);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  i = style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushstring(sty->L, gt_str_get(value));
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
}

GtStyleQueryStatus gt_style_get_num(const GtStyle *sty, const char *section,
                                    const char *key, double *val,
                                    GtFeatureNode *gn, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && key && section && val);
  gt_error_check(err);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  /* get section */
  i = style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1))
  {
    if (gn) {
      GtGenomeNode *gn_lua = gt_genome_node_ref((GtGenomeNode*) gn);
      gt_lua_genome_node_push(sty->L, gn_lua);
    }
    if (lua_pcall(sty->L, gn ? 1 : 0, 1, 0) != 0) {
      gt_error_set(err, "%s", lua_tostring(sty->L, -1));
      lua_pop(sty->L, 3);
      gt_assert(lua_gettop(sty->L) == stack_size);
      gt_rwlock_unlock(sty->lock);
      return GT_STYLE_QUERY_ERROR;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_isnumber(sty->L, -1)) {
    lua_pop(sty->L, i+1);
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  } else i++;
  /* retrieve value */
  *val = lua_tonumber(sty->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return GT_STYLE_QUERY_OK;
}

void gt_style_set_num(GtStyle *sty, const char *section, const char *key,
                    double number)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && section && key);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  i = style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushnumber(sty->L, number);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
}

GtStyleQueryStatus gt_style_get_bool(const GtStyle *sty, const char *section,
                                     const char *key, bool *val,
                                     GtFeatureNode *gn, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && key && section);
  gt_error_check(err);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  /* get section */
  i = style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1))
  {
    if (gn) {
      GtGenomeNode *gn_lua = gt_genome_node_ref((GtGenomeNode*) gn);
      gt_lua_genome_node_push(sty->L, gn_lua);
    }
    if (lua_pcall(sty->L, gn ? 1 : 0, 1, 0) != 0)
    {
      gt_error_set(err, "%s", lua_tostring(sty->L, -1));
      lua_pop(sty->L, 3);
      gt_assert(lua_gettop(sty->L) == stack_size);
      gt_rwlock_unlock(sty->lock);
      return GT_STYLE_QUERY_ERROR;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_isboolean(sty->L, -1)) {
    lua_pop(sty->L, i+1);
    gt_assert(lua_gettop(sty->L) == stack_size);
    gt_rwlock_unlock(sty->lock);
    return GT_STYLE_QUERY_NOT_SET;
  } else i++;

  /* retrieve value */
  *val = lua_toboolean(sty->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return GT_STYLE_QUERY_OK;
}

void gt_style_set_bool(GtStyle *sty, const char *section, const char *key,
                       bool val)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int i = 0;
  gt_assert(sty && section && key);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  i = style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushboolean(sty->L, val);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
}

void gt_style_unset(GtStyle *sty, const char *section, const char *key)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(sty && section && key);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  lua_getglobal(sty->L, "style");
  if (!lua_isnil(sty->L, -1)) {
    gt_assert(lua_istable(sty->L, -1));
    lua_getfield(sty->L, -1, section);
    if (!lua_isnil(sty->L, -1)) {
      gt_assert(lua_istable(sty->L, -1));
      lua_pushstring(sty->L, key);
      lua_pushnil(sty->L);
      lua_settable(sty->L, -3);
    }
    lua_pop(sty->L, 1);
  }
  lua_pop(sty->L, 1);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
}

int gt_style_to_str(const GtStyle *sty, GtStr *outstr, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int had_err;
  gt_error_check(err);
  gt_assert(sty && outstr);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);
#endif
  lua_getglobal(sty->L, "style");
  gt_str_append_cstr(outstr, "style = {\n");
  if (lua_istable(sty->L, -1))
    had_err = gt_lua_table_to_str(sty->L, outstr, -1, err);
  else {
    gt_error_set(err, "'style' must be a table. Check whether a top-level"
                      "table of this name exists.");
    had_err = -1;
  }
  gt_str_append_cstr(outstr, "}");
  lua_pop(sty->L, 1);
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return had_err;
}

int gt_style_load_str(GtStyle *sty, GtStr *instr, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  int had_err = 0;
  gt_error_check(err);
  gt_assert(sty && instr);
  gt_rwlock_wrlock(sty->lock);
#ifndef NDEBUG
  stack_size = lua_gettop(sty->L);;
#endif
  if (luaL_loadbuffer(sty->L, gt_str_get(instr), gt_str_length(instr), "str") ||
      lua_pcall(sty->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run style buffer: %s",
              lua_tostring(sty->L, -1));
    had_err = -1;
    lua_pop(sty->L, 1);
  }
  gt_assert(lua_gettop(sty->L) == stack_size);
  gt_rwlock_unlock(sty->lock);
  return had_err;
}

GtStyle* gt_style_clone(const GtStyle *sty, GtError *err)
{
  int had_err = 0;
  GtStr *sty_buffer = gt_str_new();
  GtStyle *new_sty;
  gt_assert(sty);
  if (!(new_sty = gt_style_new(err)))
    had_err = -1;
  gt_rwlock_wrlock(sty->clone_lock);
  if (!had_err)
    had_err = gt_style_to_str(sty, sty_buffer, err);
  if (!had_err)
    had_err = gt_style_load_str(new_sty, sty_buffer, err);
  gt_rwlock_unlock(sty->clone_lock);
  gt_str_delete(sty_buffer);
  return new_sty;
}

int gt_style_unit_test(GtError *err)
{
  int had_err = 0;
  GtStyle *sty = NULL, *new_sty = NULL;
  bool val = false;
  GtError *testerr;
  GtStr *test1      = gt_str_new_cstr("mRNA"),
        *str        = gt_str_new(),
        *sty_buffer = gt_str_new();
  GtColor col1, col2, col, defcol, tmpcol;
  double num = 10.0;
  gt_error_check(err);

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;col1.alpha=0.5;
  col2.red=.4;col2.green=.5;col2.blue=.6;col2.alpha=0.5;
  col.red=1.0;col.green=1.0;col.blue=1.0;col.alpha=0.5;
  defcol.red=.5;defcol.green=.5;defcol.blue=.5;defcol.alpha=0.5;

  testerr = gt_error_new();

  /* instantiate new style object */
  ensure(had_err, (sty = gt_style_new(testerr)) != NULL);
  ensure(had_err, !gt_error_is_set(testerr));

  /* test safe/unsafe mode switching */
  ensure(had_err, !gt_style_is_unsafe(sty));
  gt_style_unsafe_mode(sty);
  ensure(had_err, gt_style_is_unsafe(sty));
  gt_style_safe_mode(sty);
  ensure(had_err, !gt_style_is_unsafe(sty));

  /* at the beginning, all values are defaults, since nothing is defined */
  ensure(had_err, gt_style_get_color(sty, "exon", "fill", &tmpcol, NULL,
                                     testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, gt_style_get_color(sty, "cds", "fill", &tmpcol, NULL,
                                     testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, gt_style_get_color(sty, "foo", "fill", &tmpcol, NULL,
                                     testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, gt_style_get_num(sty, "format", "margins", &num, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, num == 10.0);
  ensure(had_err, gt_style_get_bool(sty, "exon", "collapse_to_parent", &val,
                                    NULL, testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, !val);

  /* change some values... */
  (void) gt_style_set_color(sty, "exon", "fill", &col1);
  gt_style_set_num(sty, "format", "margins", 11.0);
  gt_style_set_num(sty, "format", "foo", 2.0);

  /* is it saved correctly? */
  (void) gt_style_get_color(sty, "exon", "fill", &tmpcol, NULL, testerr);
  ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, gt_color_equals(&tmpcol, &col1));
  ensure(had_err, gt_style_get_num(sty, "format", "margins", &num, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, num == 11.0);
  ensure(had_err, !gt_error_is_set(testerr));
  ensure(had_err, gt_style_get_num(sty, "format", "foo", &num, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, num == 2.0);
  ensure(had_err, !gt_error_is_set(testerr));

  /* create a new color definition */
  gt_style_set_color(sty, "foo", "fill", &col2);
  gt_style_set_str(sty, "bar", "baz", test1);

  /* is it saved correctly? */
  ensure(had_err, gt_style_get_color(sty, "foo", "fill", &tmpcol, NULL,
                                     testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, gt_color_equals(&tmpcol, &col2));
  gt_str_reset(str);
  ensure(had_err, gt_style_get_str(sty, "bar", "baz", str, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, (strcmp(gt_str_get(str),"")!=0));
  ensure(had_err, (gt_str_cmp(str,test1)==0));
  gt_str_reset(str);
  ensure(had_err, gt_style_get_str(sty, "bar", "test", str, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, (strcmp(gt_str_get(str),"")==0));

  /* clone a GtStyle object */
  new_sty = gt_style_clone(sty, testerr);
  ensure(had_err, new_sty  != NULL);
  ensure(had_err, !gt_error_is_set(testerr));

  /* check again */
  ensure(had_err, gt_style_get_color(new_sty, "foo", "fill", &tmpcol, NULL,
                                     testerr) == 0);
  ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, gt_color_equals(&tmpcol, &col2));
  gt_str_reset(str);
  ensure(had_err, gt_style_get_str(new_sty, "bar", "baz", str, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, (strcmp(gt_str_get(str),"")!=0));
  ensure(had_err, (gt_str_cmp(str,test1)==0));
  gt_str_reset(str);
  ensure(had_err, gt_style_get_str(new_sty, "bar", "test", str, NULL,
                                   testerr) != GT_STYLE_QUERY_ERROR);
  ensure(had_err, (strcmp(gt_str_get(str),"")==0));

  /* mem cleanup */
  gt_error_delete(testerr);
  gt_str_delete(test1);
  gt_str_delete(str);
  gt_str_delete(sty_buffer);
  gt_style_delete(sty);
  gt_style_delete(new_sty);

  return had_err;
}

void gt_style_delete_without_state(GtStyle *sty)
{
  if (!sty) return;
  gt_rwlock_wrlock(sty->lock);
  if (sty->reference_count)
  {
    sty->reference_count--;
    gt_rwlock_unlock(sty->lock);
    return;
  }
  gt_free(sty->filename);
  gt_rwlock_unlock(sty->lock);
  gt_rwlock_delete(sty->lock);
  gt_rwlock_delete(sty->clone_lock);
  gt_free(sty);
}

void gt_style_delete(GtStyle *style)
{
  if (!style) return;
  gt_rwlock_wrlock(style->lock);
  if (style->reference_count)
  {
    style->reference_count--;
    gt_rwlock_unlock(style->lock);
    return;
  }
  if (style->L) lua_close(style->L);
  gt_rwlock_unlock(style->lock);
  gt_style_delete_without_state(style);
}
