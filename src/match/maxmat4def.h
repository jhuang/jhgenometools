/*
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef MAXMAT4DEF_H
#define MAXMAT4DEF_H

#include "match/matchmode_api.h"
#include "core/array_api.h"
#include "core/alphabet_api.h"
#include "core/encseq_api.h"
#include "core/defined-types.h"

/*
  This file defines some constants and types for computing maximal 
  matches using fm-index.
*/

/*
  The following structure contains all information 
  specified for the function showmatch
*/
typedef struct
{
	bool showstring,
       showreversepositions,
       showsequencelengths;
  enum GtReadmode queryreadmode;
  GtArray *mumcandtab;
} Showspecinfo;

/*
  The following structure contains all information
  required while computing and processing the matches.
*/
typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  GtMatchmode matchmode;
  const GtAlphabet *alphabet;
  //Preprocessmatchfunction preprocessmatchfunction;
  //Processmatchfunction processmatchfunction;
  //Postprocessmatchfunction postprocessmatchfunction;
  //Showmatchfunction showmatchfunction;  
  const GtEncseq *encseq;
  Definedunsignedlong leastlength;
  Showspecinfo *showspecinfo;
} Matchprocessinfo;

/*
  Functions processing a maximal match are of the following type.
*/

typedef short int (*Processmatchfunction)
                            (const GtEncseq *encseq,
                             const GtUchar *query,
                             unsigned long querypos,
                             unsigned long querylength,
                             unsigned long matchlength,
                             unsigned long subjectpos,
                             Showspecinfo *showspecinfo); 

#endif
