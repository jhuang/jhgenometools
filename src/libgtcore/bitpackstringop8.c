/*
** autogenerated content - DO NOT EDIT
*/
/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**
** See LICENSE file or http://genometools.org/license.html for license details.
**
*/
#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>

#include "bitpackstring.h"

#define MIN(a, b) (((a)<(b))?(a):(b))
#define MIN3(a, b, c) (((a)<(b))?((a)<(c)?(a):(c)):((b)<(c)?(b):(c)))

uint8_t
bsGetUInt8(const BitString str, BitOffset offset, unsigned numBits)
{
  uint_fast32_t accum = 0;
  unsigned bitsLeft = numBits, bitTop = offset%bitElemBits;
  size_t elemStart = offset/bitElemBits;
  const BitElem *p = str + elemStart;
  assert(str);
#ifndef NDEBUG
  if(numBits > sizeof(accum)*CHAR_BIT)
    fprintf(stderr, "numBits = %u\n", numBits);
#endif
  assert(numBits <= sizeof(accum)*CHAR_BIT);
  if(bitTop)
  {
    uint_fast32_t mask;
    unsigned bits2Read = MIN(bitElemBits - bitTop, bitsLeft);
    unsigned unreadRightBits = (bitElemBits - bitTop - bits2Read);
    mask = (~((~(uint_fast32_t)0) << bits2Read)) << unreadRightBits;
    accum = ((*p++) & mask) >> unreadRightBits;
    bitsLeft -= bits2Read;
  }
  /* get bits from intervening elems */
  while(bitsLeft >= bitElemBits)
  {
    accum = accum << bitElemBits | (*p++);
    bitsLeft -= bitElemBits;
  }
  /* get bits from last elem */
  accum = accum << bitsLeft |
    (((*p) & ((~(uint_fast32_t)0)<<(bitElemBits - bitsLeft)))
     >>(bitElemBits - bitsLeft));
  return accum;
}

void
bsStoreUInt8(BitString str, BitOffset offset,
                 unsigned numBits, uint8_t val)
{
  unsigned bitsLeft = numBits,
    bitTop = offset%bitElemBits;
  size_t elemStart = offset/bitElemBits;
  BitElem *p = str + elemStart;
  assert(str);
  assert(numBits <= sizeof(val)*CHAR_BIT);
  /* set bits of first element, accounting for bits to be preserved */
  if(bitTop)
  {
    uint_fast32_t mask = ~(uint_fast32_t)0;
    if(bitElemBits < (sizeof(uint_fast32_t)*CHAR_BIT))
    {
      mask <<= bitElemBits;
    }
    else
    {
      mask = 0;
    }
    mask = (~mask) >> bitTop;
    if(numBits < bitElemBits - bitTop)
    {
      unsigned backShift = bitElemBits - numBits - bitTop;
      mask &= ~(uint_fast32_t)0 << backShift;
      *p = (*p & ~mask) | ((val << backShift) & mask);
      /* TODO: try wether  r = a ^ ((a ^ b) & mask) is faster, see below */
      return;
    }
    else
    {
      bitsLeft -= bitElemBits - bitTop;
      *p = (*p & ~mask) | ((val >> bitsLeft) & mask);
      ++p;
    }
  }
  /* set bits for intervening elems */
  while(bitsLeft >= bitElemBits)
  {
    bitsLeft -= bitElemBits;
    *p++ = val >> bitsLeft;
  }
  /* set bits for last elem */
  {
    uint_fast32_t mask = ((~(uint_fast32_t)0)<<(bitElemBits - bitsLeft));
    if(bitElemBits < (sizeof(uint8_t)*CHAR_BIT))
      mask &= (~(~(uint_fast32_t)0<<bitElemBits));
    *p = (*p & ~mask) | ((val << (bitElemBits - bitsLeft)) & mask);
  }
}

/**************************************************************************/
/* Merge bits from two values according to a mask                         */
/*                                                                        */
/* unsigned int a;    // value to merge in non-masked bits                */
/* unsigned int b;    // value to merge in masked bits                    */
/* unsigned int mask; // 1 where bits from b should be selected;          */
/*                    // 0 where from a.                                  */
/* unsigned int r;    // result of (a & ~mask) | (b & mask) goes here     */
/*                                                                        */
/* r = a ^ ((a ^ b) & mask);                                              */
/*                                                                        */
/* This shaves one operation from the obvious way of combining two sets   */
/* of bits according to a bit mask. If the mask is a constant, then there */
/* may be no advantage.                                                   */
/*                                                                        */
/* Ron Jeffery sent this to me on February 9, 2006.                       */
/**************************************************************************/

void
bsGetUniformUInt8Array(const BitString str, BitOffset offset,
                           unsigned numBits, size_t numValues,
                           uint8_t val[])
{
  /* idea: read as much as possible from str in each iteration,
   * accumulate if bitsLeft < numBits */
  size_t j = 0;
  BitOffset totalBitsLeft = numValues * numBits;
  size_t elemStart = offset/bitElemBits;
  unsigned bitTop = offset%bitElemBits,
    bitsRead = 0; /*< how many bits in current *p are read */
  const BitElem *p = str + elemStart;
  unsigned bitsInAccum = 0;
  uint_fast32_t accum = 0, valMask = ~(uint_fast32_t)0;
  if(numBits < (sizeof(val[0])*CHAR_BIT))
    valMask = ~(valMask << numBits);
  assert(str && val);
  assert(numBits <= sizeof(val[0])*CHAR_BIT);
  /* user requested zero values, ugly but must be handled, since legal */
  if(!totalBitsLeft)
  {
    return;
  }
  /* get bits of first element if not aligned */
  if(bitTop)
  {
    uint_fast32_t mask; /*< all of the bits we want to get from *p */
    unsigned bits2Read = MIN(bitElemBits - bitTop, totalBitsLeft);
    unsigned unreadRightBits = (bitElemBits - bitTop - bits2Read);
    mask = (~((~(uint_fast32_t)0) << bits2Read)) << unreadRightBits;
    accum = ((*p++) & mask) >> unreadRightBits;
    bitsInAccum += bits2Read;
    totalBitsLeft -= bits2Read;
  }
  while(j < numValues)
  {
    while(bitsInAccum < numBits && totalBitsLeft)
    {
      unsigned bits2Read, bitsFree = sizeof(accum)*CHAR_BIT - bitsInAccum;
      uint_fast32_t mask;
      bits2Read = MIN3(bitsFree, bitElemBits - bitsRead, totalBitsLeft);
      mask = (~((~(uint_fast32_t)0) << bits2Read));
      accum = accum << bits2Read | (((*p) >> (bitElemBits
                                              - bits2Read - bitsRead)) & mask);
      bitsInAccum += bits2Read;
      totalBitsLeft -= bits2Read;
      /* all of *p consumed? */
      if((bitsRead += bits2Read) == bitElemBits)
      {
        ++p, bitsRead = 0;
      }
    }
    /* now we have enough bits in accum */
    while(bitsInAccum >= numBits)
    {
      val[j++] = ((accum >> (bitsInAccum - numBits)) & valMask );
      bitsInAccum -= numBits;
    }
  }
}

void
bsStoreUniformUInt8Array(BitString str, BitOffset offset, unsigned numBits,
                             size_t numValues, const uint8_t val[])
{
  /* idea: read as much as possible from val in each iteration,
   * accumulate if bitsInAccum < bitElemBits */
  size_t j = 0;
  BitOffset totalBitsLeft = numValues * numBits;
  unsigned bitTop = offset%bitElemBits,
    bitsLeft; /*< how many bits in currentVal == val[j] are left */
  BitElem *p = str + offset/bitElemBits;
  unsigned bitsInAccum;
  uint_fast32_t accum, valMask = ~(uint_fast32_t)0, currentVal;
  if(numBits < (sizeof(val[0])*CHAR_BIT))
    valMask = ~(valMask << numBits);
  assert(str && val);
  assert(numBits <= sizeof(val[0])*CHAR_BIT);
  /* user requested zero values, ugly but must be handled, since legal */
  if(!totalBitsLeft)
  {
    return;
  }
  accum = val[0] & valMask;
  totalBitsLeft -= bitsInAccum = numBits;
  if(totalBitsLeft)
  {
    currentVal = val[++j] & valMask;
    totalBitsLeft -= bitsLeft = numBits;
  }
  else
  {
    currentVal = 0;
    bitsLeft = 0;
  }
  /* set bits of first element if not aligned */
  if(bitTop)
  {
    BitElem mask = ~(~(uint_fast32_t)0 << (bitElemBits - bitTop));
    while((totalBitsLeft || bitsLeft) && bitsInAccum < bitElemBits - bitTop)
    {
      unsigned bits2Read, bitsFree = sizeof(accum)*CHAR_BIT - bitsInAccum;

      if((bits2Read = MIN(bitsFree, bitsLeft)) < sizeof(accum)*CHAR_BIT)
      {
        accum = accum << bits2Read
          | ((currentVal) >> (bitsLeft - bits2Read));
      }
      else
        accum = currentVal;

      /* all of val[j] consumed? */
      bitsInAccum += bits2Read;
      if(!(bits2Read -= bitsLeft))
        currentVal = val[++j] & valMask, totalBitsLeft -= bitsLeft = numBits;
    }
    /* at this point accum holds as many bits as we could get
     * to fill the first BitElem in str, but did we get enough? */
    if(bitsInAccum < bitElemBits - bitTop)
    {
      /* no there's not enough */
      unsigned backShift = bitElemBits - bitsInAccum - bitTop;
      mask &= ~(uint_fast32_t)0 << backShift;
      *p = (*p & ~mask) | ((accum << backShift) & mask);
      /* TODO: try wether  r = a ^ ((a ^ b) & mask) is faster, see below */
      return; /* if we couldn't gather more bits, there's none left */
    }
    else
    {
      /* yep, just or with accumVals */
      *p = (*p & ~mask) | (accum >> (bitsInAccum - bitElemBits + bitTop));
      ++p;
      bitsInAccum -= bitElemBits - bitTop;
    }
  }

  while(totalBitsLeft || (bitsInAccum + bitsLeft) > bitElemBits)
  {
    while((totalBitsLeft || bitsLeft)
          && ((bitsInAccum < bitElemBits)
              || (bitsLeft < sizeof(accum)*CHAR_BIT - bitsInAccum)))
    {
      unsigned bits2Read, bitsFree = sizeof(accum)*CHAR_BIT - bitsInAccum;
      if((bits2Read = MIN(bitsFree, bitsLeft)) < sizeof(accum)*CHAR_BIT)
      {
        uint_fast32_t mask = ~((~(uint_fast32_t)0) << bits2Read);
        accum = accum << bits2Read
          | ((currentVal >> (bitsLeft - bits2Read)) & mask);
      }
      else
        accum = currentVal;
      bitsInAccum += bits2Read;
      /* all of currentVal == val[j] consumed? */
      if(bits2Read == bitsLeft && totalBitsLeft)
        currentVal = val[++j] & valMask, totalBitsLeft -= bitsLeft = numBits;
      else
        bitsLeft -= bits2Read;
    }
    /* now we have enough bits in accum */
    while(bitsInAccum >= bitElemBits)
    {
      *p++ = accum >> (bitsInAccum - bitElemBits);
      bitsInAccum -= bitElemBits;
    }
  }
  /* write the rest bits left in accum and currentVal */
  accum = (accum << bitsLeft)
    | (currentVal & (valMask >> (numBits - bitsLeft)));
  bitsInAccum += bitsLeft;
  while(bitsInAccum >= bitElemBits)
  {
    *p++ = accum >> (bitsInAccum - bitElemBits);
    bitsInAccum -= bitElemBits;
  }
  if(bitsInAccum)
  {
    uint_fast32_t mask = ~(uint_fast32_t)0 << (bitElemBits - bitsInAccum);
    *p = (*p & ~mask) | ((accum << (bitElemBits - bitsInAccum))& mask);
  }
}
