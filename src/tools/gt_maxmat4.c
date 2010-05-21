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


//#include "core/defined-types.h"
//#include "core/error.h"
//#include "core/option.h"
//#include "core/unused_api.h"
//#include "core/versionfunc.h"
//#include "match/sarr-def.h"
#include "match/maxmat4.h"
#include "match/esa-map.h"
#include "match/eis-voiditf.h"
#include "tools/gt_maxmat4.h"



typedef struct {
  bool mum,
       mumcand,
       nucleotidesonly,                      // matchnucleotidesonly
       bothdirections,                       // computebothdirections
       reversecomplement,                      // onlyreversecomplementmatches
       showstring,                             // showstring
       showreversepositions,                     // showreversepositions
       showsequencelengths;                      // showsequencelengths
  Definedunsignedlong leastlength;               // leastlength  
  bool verbose;
} GtMaxmat4Arguments;



static void* gt_maxmat4_arguments_new(void)
{
  GtMaxmat4Arguments *arguments;
  arguments = gt_malloc(sizeof (GtMaxmat4Arguments));
  return arguments;
}

static void gt_maxmat4_arguments_delete(void *tool_arguments)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  if (!arguments)
  {
    return;
  }
  gt_free(arguments);
}

static GtOptionParser* gt_maxmat4_option_parser_new(void *tool_arguments)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  GtOption *option_mum, *option_mumcand, *option_n, *option_l, 
         *option_b, *option_r, *option_s, *option_c, *option_L;
  GtOptionParser *op;

  gt_assert(arguments);
  
  op = gt_option_parser_new("[options] <reference-file> <query-files>[...]",
                            "Find and output (to stdout) the positions and length of all "
                            "sufficiently long (unique) maximal matches of a substring in "
                            "<reference-file>(in format of specify packed index) and <query-file>.");
  gt_option_parser_set_mailaddress(op,"<gt-users@genometools.org>");

  /* -mum */
  option_mum = gt_option_new_bool("mum", 
                                  "compute MUMs, i.e. maximal matches that are unique "
                                  "in both sequences ",
                                  &arguments->mum, false);
  gt_option_parser_add_option(op, option_mum);

  /* -mumcand */
  option_mumcand = gt_option_new_bool("mumcand", 
                                      "compute MUM-candidates, i.e. maximal matches that are unique "
                                      "in the subject-sequence but not necessarily in the query-sequence ", 
                                      &arguments->mumcand, false);
  gt_option_parser_add_option(op, option_mumcand);

  /* -n for matchnucleotidesonly */
  option_n = gt_option_new_bool("n", 
                                "match only the characters a, c, g, or t "
                                "they can be in upper or in lower case",
                                &arguments->nucleotidesonly, false);
  gt_option_parser_add_option(op, option_n);
  
  /* -l for leastlength */
  option_l = gt_option_new_ulong_min("l",
                                 "set the minimum length of a match "
                                 "if not set, the default value is 20",
                                 &arguments->leastlength.valueunsignedlong,
                                 2,(unsigned long) 1);
  gt_option_parser_add_option(op, option_l);
  
  /* -b for computebothdirections */
  option_b = gt_option_new_bool("b",
                                "compute forward and reverse complement matches", 
                                &arguments->bothdirections, false);
  gt_option_parser_add_option(op, option_b);

  /* -r for onlyreversecomplementmatches */
  option_r = gt_option_new_bool("r", 
                                "only compute reverse complement matches",
                                &arguments->reversecomplement, false);
  gt_option_parser_add_option(op, option_r);

  /* -s for showstring */
  option_s = gt_option_new_bool("s",
                                "show the matching substrings",
                                &arguments->showstring, false);
  gt_option_parser_add_option(op, option_s);

  /* -c for showreversepositions */
  option_c = gt_option_new_bool("c", 
                                "report the query-position of a reverse complement match "
                                "relative to the original query sequence",
                                &arguments->showreversepositions, false);
  gt_option_parser_add_option(op, option_c);

  /* -L for showsequencelengths */
  option_L = gt_option_new_bool("L", 
                                "show the length of the query sequences on the header line",
                                &arguments->showsequencelengths, false);
  gt_option_parser_add_option(op, option_L);


  /* option implications */
  gt_option_imply_either_2(option_c, option_b, option_r);

  /* option exclusions */
  gt_option_exclude(option_mum, option_mumcand);
  gt_option_exclude(option_b, option_r);

  /* set minimal arugments */
  gt_option_parser_set_min_args(op, 2);

  return op;
}

static int gt_maxmat4_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments,
                                      GtError *err)
{
  GT_UNUSED GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  int had_err = 0;
  
  gt_error_check(err);
  gt_assert(arguments);


  return had_err;
}

static int gt_maxmat4_runner(GT_UNUSED int argc,
                             GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args,
                             void *tool_arguments, 
                             GT_UNUSED GtError *err)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  //GtMaxmat4 *maxmat4;
  int arg = parsed_args;

  Suffixarray suffixarray;
  void *packedindex = NULL;
  GtLogger *logger = NULL;
  bool haserr = false;
  const GtAlphabet *alphabet = NULL;
  // unsigned int prefixlength = 0;
  unsigned long totallength;
  
  // init the referencefile and queryfiles
  GtStr *referencefile = gt_str_new_cstr(argv[arg]);
  /* it makes only sense if we got a reference file in format pck */
  //if (gt_option_is_set(optionpckindex))
  //{
      //gfmsubcallinfo->indextype = Packedindextype;
  //} else
  //{
	//gt_error_set(err,"a file must in format of specify packed index for <reference file> be used");
	//had_err = -1;
  //}
  GtStrArray *queryfiles = gt_str_array_new();
  unsigned short idx;
  for (idx=arg+1; idx<argc; idx++)
  {
	gt_str_array_add_cstr(queryfiles, argv[idx]); 
  }
  
  /* load packed index in memory */
  gt_error_check(err);
  gt_assert(tool_arguments);
  logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
  unsigned int mappedbits = SARR_ESQTAB;
  
  // map suffixarray from referencefile, a referencefile contains maybe many sequence, is this merged in a suffixarray?
	if (gt_mapsuffixarray(&suffixarray,
						 mappedbits,
						 referencefile,
						 logger,
						 err) != 0)
	{
		haserr = true;
		totallength = 0;
	} else
	{
		alphabet = gt_encseq_alphabet(suffixarray.encseq);
		// prefixlength = suffixarray.prefixlength;
		totallength = gt_encseq_total_length(suffixarray.encseq);
	}

	if (!haserr)
	{
		// load packed index from suffxiarray
		packedindex = gt_loadvoidBWTSeqForSA(referencefile,
											&suffixarray,
											totallength,
											false,
											err);
		if (packedindex == NULL)
		{
			haserr = true;
		}
	}


	if (!haserr)
	{
		int readmode = GT_READMODE_FORWARD;  // default value is forward
		if (arguments->reversecomplement) 
		{
			readmode = GT_READMODE_REVCOMPL;
			//printf ("%s \n", "reversecomplement is true, correspondent method is still in work");
		}
		else if (arguments->bothdirections) 
		{
			readmode = 4;  // for both directions
			//printf ("%s \n", "bothdirections is true, correspondent method is still in work");
		}
		else  
		{
			if (gt_findmum(suffixarray.encseq,
										packedindex,
										totallength,
										alphabet,
										queryfiles,
										readmode,
										arguments->leastlength,
										arguments->nucleotidesonly,
										//arguments->bothdirections,
										//arguments->reversecomplement,
										arguments->showstring,
										arguments->showreversepositions,
										arguments->showsequencelengths,
										arguments->verbose,
										err) != 0)
			{
				haserr = true;
			}
		}
	}
 

  /* maxmat4 construction */
  //maxmat4 = gt_maxmat4_new(referencefile, queryfiles, err);
  //if (!maxmat4)
    //had_err = -1;

  /* free */
  //gt_maxmat4_delete(maxmat4);

  if (packedindex != NULL)
  {
    gt_deletevoidBWTSeq(packedindex);
  }
  gt_freesuffixarray(&suffixarray);
  
  gt_logger_delete(logger);
  gt_str_delete(referencefile);
  gt_str_array_delete(queryfiles);
  return haserr ? -1 : 0;
}

GtTool* gt_maxmat4(void)
{
  return gt_tool_new(gt_maxmat4_arguments_new,
                     gt_maxmat4_arguments_delete,
                     gt_maxmat4_option_parser_new,
                     gt_maxmat4_arguments_check,
                     gt_maxmat4_runner);
}
