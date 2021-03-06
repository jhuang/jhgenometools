/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/types_api.h"
#include "core/unused_api.h"

#include "spacedef.h"
#include "esa-seqread.h"
#include "esa-dfs.h"
#include "stamp.h"

#define ABOVETOP  stackspace[nextfreeItvinfo]
#define TOP       stackspace[nextfreeItvinfo-1]
#define BELOWTOP  stackspace[nextfreeItvinfo-2]

#define INCSTACKSIZE  32

#define PUSHDFS(D,B,PREVIOUSPTR)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocItvinfo(PREVIOUSPTR,\
                                    allocatedItvinfo,\
                                    allocatedItvinfo+INCSTACKSIZE,\
                                    allocateDfsinfo,\
                                    state);\
          allocatedItvinfo += INCSTACKSIZE;\
        }\
        stackspace[nextfreeItvinfo].depth = D;\
        stackspace[nextfreeItvinfo].lastisleafedge = B;\
        nextfreeItvinfo++

typedef struct
{
  bool lastisleafedge;
  unsigned long depth;
  Dfsinfo *dfsinfo;
} Itvinfo;

static Itvinfo *allocItvinfo(Itvinfo *ptr,
                             unsigned long currentallocated,
                             unsigned long allocated,
                             Dfsinfo *(*allocateDfsinfo)(Dfsstate *),
                             Dfsstate *state)
{
  unsigned long i;
  Itvinfo *itvinfo;

  ALLOCASSIGNSPACE(itvinfo,ptr,Itvinfo,allocated);
  if (allocateDfsinfo != NULL)
  {
    gt_assert(allocated > currentallocated);
    for (i=currentallocated; i<allocated; i++)
    {
      itvinfo[i].dfsinfo = allocateDfsinfo(state);
    }
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

static void freeItvinfo(Itvinfo *ptr,
                        unsigned long allocated,
                        void (*freeDfsinfo)(Dfsinfo *,Dfsstate *),
                        Dfsstate *state)
{
  unsigned long i;

  for (i=0; i<allocated; i++)
  {
    freeDfsinfo(ptr[i].dfsinfo,state);
  }
  FREESPACE(ptr);
}

int gt_depthfirstesa(Sequentialsuffixarrayreader *ssar,
                  Dfsinfo *(*allocateDfsinfo)(Dfsstate *),
                  void(*freeDfsinfo)(Dfsinfo *,Dfsstate *),
                  int(*processleafedge)(bool,unsigned long,Dfsinfo *,
                                        unsigned long,Dfsstate *,
                                        GtError *),
                  int(*processbranchedge)(bool,
                                          unsigned long,
                                          Dfsinfo *,
                                          Dfsinfo *,
                                          Dfsstate *,
                                          GtError *),
                  int(*processcompletenode)(unsigned long,
                                            Dfsinfo *,unsigned long,
                                            Dfsstate *,GtError *),
                  void (*assignleftmostleaf)(Dfsinfo *,unsigned long,
                                             Dfsstate *),
                  void (*assignrightmostleaf)(Dfsinfo *,unsigned long,
                                              unsigned long,
                                              unsigned long,Dfsstate *),
                  Dfsstate *state,
                  GT_UNUSED GtLogger *logger,
                  GtError *err)
{
  bool firstedge,
       firstrootedge;
  unsigned long previoussuffix = 0,
         previouslcp,
         currentindex,
         currentlcp = 0; /* May be necessary if currentlcp is used after the
                            outer while loop */
  unsigned long allocatedItvinfo = 0,
                nextfreeItvinfo = 0;
  Itvinfo *stackspace;
  bool haserr = false;

#ifdef INLINEDSequentialsuffixarrayreader
  GtUchar tmpsmalllcpvalue;
  gt_logger_log(logger,"# inlined Sequentialsuffixarrayreader\n");
#else
  int retval;
#endif
  firstrootedge = true;
  PUSHDFS(0,true,NULL);
  if (assignleftmostleaf != NULL)
  {
    assignleftmostleaf(TOP.dfsinfo,0,state);
  }
  for (currentindex = 0; !haserr; currentindex++)
  {
#ifdef INLINEDSequentialsuffixarrayreader
    NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
#else
    retval = gt_nextSequentiallcpvalue(&currentlcp,ssar,err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    retval = gt_nextSequentialsuftabvalue(&previoussuffix,ssar);
    gt_assert(retval >= 0);
    if (retval == 0)
    {
      gt_error_set(err,"Missing value in suftab");
      haserr = true;
      break;
    }
#endif
    while (currentlcp < TOP.depth)
    {
      if (TOP.lastisleafedge)
      {
        if (processleafedge != NULL &&
            processleafedge(false,TOP.depth,TOP.dfsinfo,
                            previoussuffix,state,err) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        gt_assert(nextfreeItvinfo < allocatedItvinfo);
        gt_assert(ABOVETOP.dfsinfo != NULL);
        if (processbranchedge != NULL &&
            processbranchedge(false,
                              TOP.depth,
                              TOP.dfsinfo,
                              ABOVETOP.dfsinfo,
                              state,
                              err) != 0)
        {
          haserr = true;
          break;
        }
      }
      if (assignrightmostleaf != NULL)
      {
        assignrightmostleaf(TOP.dfsinfo,
                            currentindex,
                            previoussuffix,
                            currentlcp,
                            state);
      }
      if (processcompletenode != NULL &&
          processcompletenode(TOP.depth,TOP.dfsinfo,
                              nextfreeItvinfo >= 2UL ? BELOWTOP.depth
                                                     : 0,
                              state,err) != 0)
      {
        haserr = true;
        break;
      }
      gt_assert(nextfreeItvinfo > 0);
      nextfreeItvinfo--;
    }
    if (haserr)
    {
      break;
    }
    if (currentlcp == TOP.depth)
    {
      if (firstrootedge && TOP.depth == 0)
      {
        firstedge = true;
        firstrootedge = false;
      } else
      {
        firstedge = false;
      }
      if (TOP.lastisleafedge)
      {
        if (processleafedge != NULL &&
            processleafedge(firstedge,TOP.depth,TOP.dfsinfo,
                            previoussuffix,state,err) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if (!firstedge)
        {
          gt_assert(nextfreeItvinfo < allocatedItvinfo);
        }
        if (processbranchedge != NULL &&
            processbranchedge(firstedge,
                              TOP.depth,
                              TOP.dfsinfo,
                              /*firstedge ? NULL : ABOVETOP.dfsinfo,*/
                              ABOVETOP.dfsinfo,
                              state,
                              err) != 0)
        {
          haserr = true;
          break;
        }
        TOP.lastisleafedge = true;
      }
    } else
    {
      PUSHDFS(currentlcp,true,stackspace);
      if (BELOWTOP.lastisleafedge)
      {
       if (assignleftmostleaf != NULL)
        {
          assignleftmostleaf(TOP.dfsinfo,currentindex,state);
        }
        if (processleafedge != NULL &&
            processleafedge(true,
                            TOP.depth,
                            TOP.dfsinfo,
                            previoussuffix,
                            state,
                            err) != 0)
        {
          haserr = true;
          break;
        }
        BELOWTOP.lastisleafedge = false;
      } else
      {
        previouslcp = TOP.depth;
        if (processbranchedge != NULL &&
            processbranchedge(true,
                              previouslcp,
                              TOP.dfsinfo,
                              NULL, /* not used since firstsucc = true */
                              state,
                              err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  if (!haserr && TOP.lastisleafedge)
  {
#ifdef INLINEDSequentialsuffixarrayreader
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
#else
    retval = gt_nextSequentialsuftabvalue(&previoussuffix,ssar);
    gt_assert(retval >= 0);
    if (retval == 0)
    {
      gt_error_set(err,"Missing value in suftab");
      haserr = true;
    }
#endif
    if (!haserr)
    {
      if (processleafedge != NULL &&
          processleafedge(false,
                          TOP.depth,
                          TOP.dfsinfo,
                          previoussuffix,
                          state,
                          err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (assignrightmostleaf != NULL)
      {
        assignrightmostleaf(TOP.dfsinfo,
                            currentindex,
                            previoussuffix,
                            currentlcp,
                            state);
      }
      if (processcompletenode != NULL &&
          processcompletenode(TOP.depth,TOP.dfsinfo,
                              nextfreeItvinfo >= 2UL ? BELOWTOP.depth
                                                     : 0,
                              state,err) != 0)
      {
        haserr = true;
      }
    }
  }
  freeItvinfo(stackspace,
              allocatedItvinfo,
              freeDfsinfo,
              state);
  return haserr ? -1 : 0;
}
