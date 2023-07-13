#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "mathsupport.h"
#include "gt-alloc.h"
#include "qgram_limits.h"
#include "array_api.h"
#include "multisets.h"

/* BEGIN{created by genmultisets.py -w -m 7 -k 20}
   DO NOT EDIT */
const unsigned int multiset_weights_table_alphasize = 20;
const unsigned int multiset_weights_table_width = 7;
const unsigned int multiset_weights_table[] = {
  /* pos=0 */
  0,
  134596,
  235543,
  310156,
  364420,
  403180,
  430312,
  448876,
  461252,
  469260,
  474265,
  477268,
  478984,
  479908,
  480370,
  480580,
  480664,
  480692,
  480699,
  480700,
  /* pos=1 */
  0,
  33649,
  59983,
  80332,
  95836,
  107464,
  116032,
  122220,
  126588,
  129591,
  131593,
  132880,
  133672,
  134134,
  134386,
  134512,
  134568,
  134589,
  134595,
  134596,
  /* pos=2 */
  0,
  7315,
  13300,
  18145,
  22021,
  25081,
  27461,
  29281,
  30646,
  31647,
  32362,
  32857,
  33187,
  33397,
  33523,
  33593,
  33628,
  33643,
  33648,
  33649,
  /* pos=3 */
  0,
  1330,
  2470,
  3439,
  4255,
  4935,
  5495,
  5950,
  6314,
  6600,
  6820,
  6985,
  7105,
  7189,
  7245,
  7280,
  7300,
  7310,
  7314,
  7315,
  /* pos=4 */
  0,
  190,
  361,
  514,
  650,
  770,
  875,
  966,
  1044,
  1110,
  1165,
  1210,
  1246,
  1274,
  1295,
  1310,
  1320,
  1326,
  1329,
  1330,
  /* pos=5 */
  0,
  19,
  37,
  54,
  70,
  85,
  99,
  112,
  124,
  135,
  145,
  154,
  162,
  169,
  175,
  180,
  184,
  187,
  189,
  190,
  /* pos=6 */
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19};
/* END{created by genmultisets.py -w -m 7 -k 20} */

static size_t binom_evaluate(size_t n,size_t k)
{
  size_t j, result;

  if (k == 0 || n == k)
    return (size_t) 1;
  if (n == 0)
    return 0;
  for (result = 1UL, j=1; j<=k; j++)
  {
    result = (result * (n - k + j))/j;
  }
  return result;
}

size_t number_of_multisets(size_t alphasize,size_t qgram_length)
{
  return binom_evaluate(alphasize + qgram_length - 1,qgram_length);
}

struct MultisetsInfo
{
  unsigned int max_qgram_length,
               **weights;
};

MultisetsInfo *multisets_info_new(unsigned long alphasize)
{
  MultisetsInfo *ms_info = gt_malloc(sizeof *ms_info);
  unsigned int idx;

  assert(alphasize > 0);
  if (alphasize == multiset_weights_table_alphasize)
  {
    assert((sizeof multiset_weights_table/sizeof multiset_weights_table[0]) ==
           alphasize * multiset_weights_table_width);
    ms_info->max_qgram_length = multiset_weights_table_width;
    ms_info->weights = gt_malloc(sizeof *ms_info->weights *
                                 ms_info->max_qgram_length);
    for (idx = 0; idx < ms_info->max_qgram_length; idx++)
    {
      ms_info->weights[idx] = (unsigned int *)
                              (multiset_weights_table + idx * alphasize);
    }
  } else
  {
    unsigned int previous_sum = 0,
                 *sums = gt_malloc(sizeof *sums * alphasize);

    ms_info->max_qgram_length = 2 + 1;
    ms_info->weights = gt_malloc(sizeof *ms_info->weights *
                                 ms_info->max_qgram_length);

    ms_info->weights[0] = gt_malloc(sizeof **ms_info->weights * alphasize);
    ms_info->weights[1] = gt_malloc(sizeof **ms_info->weights * alphasize);
    for (idx = 1; idx < alphasize; idx++)
    {
      previous_sum += idx;
      sums[alphasize - idx] = previous_sum;
    }
    sums[0] = 0;
    ms_info->weights[0][0] = 0;
    for (idx = 1; idx < alphasize; idx++)
    {
      ms_info->weights[0][idx] = ms_info->weights[0][idx-1] + sums[idx];
    }
    for (idx = 0; idx < alphasize; idx++)
    {
      ms_info->weights[1][idx] = (idx * (2 * alphasize - idx - 1))/2U;
    }
    for (idx = 0; idx < alphasize; idx++)
    {
      printf("%u %u %u\n",idx,ms_info->weights[0][idx],
                              ms_info->weights[1][idx]);
    }
  }
  return ms_info;
}

unsigned int **multisets_weights(const MultisetsInfo *ms_info,
                                 unsigned long qgram_length)
{
  assert(qgram_length >= 2 && qgram_length <= MAX_QGRAM_LENGTH);
  return ms_info->weights + multiset_weights_table_width - qgram_length;
}

unsigned long multisets_encode(const MultisetsInfo *ms_info,const uint8_t *ms,
                               unsigned long qgram_length)
{
  unsigned int **weights = multisets_weights(ms_info,qgram_length);

  switch (qgram_length)
  {
    case 2:
      return weights[0][(unsigned long) ms[0]] +
             (unsigned int) ms[1];
    case 3:
      return weights[0][(unsigned long) ms[0]] +
             weights[1][(unsigned long) ms[1]] +
             (unsigned int) ms[2];
    case 4:
      return weights[0][(unsigned long) ms[0]] +
             weights[1][(unsigned long) ms[1]] +
             weights[2][(unsigned long) ms[2]] +
             (unsigned int) ms[3];
    case 5:
      return weights[0][(unsigned long) ms[0]] +
             weights[1][(unsigned long) ms[1]] +
             weights[2][(unsigned long) ms[2]] +
             weights[3][(unsigned long) ms[3]] +
             (unsigned int) ms[4];
    case 6:
      return weights[0][(unsigned long) ms[0]] +
             weights[1][(unsigned long) ms[1]] +
             weights[2][(unsigned long) ms[2]] +
             weights[3][(unsigned long) ms[3]] +
             weights[4][(unsigned long) ms[4]] +
             (unsigned int) ms[5];
    case 7:
      return weights[0][(unsigned long) ms[0]] +
             weights[1][(unsigned long) ms[1]] +
             weights[2][(unsigned long) ms[2]] +
             weights[3][(unsigned long) ms[3]] +
             weights[4][(unsigned long) ms[4]] +
             weights[5][(unsigned long) ms[5]] +
             (unsigned int) ms[6];
     default:
       fprintf(stderr,"Illegal qgram_length %lu\n",qgram_length);
       exit(EXIT_FAILURE);
  }
}

void multisets_info_delete(MultisetsInfo *ms_info)
{
  if (ms_info != NULL)
  {
    if (ms_info->max_qgram_length != multiset_weights_table_width)
    {
      gt_free(ms_info->weights[0]);
      gt_free(ms_info->weights[1]);
    }
    gt_free(ms_info->weights);
    gt_free(ms_info);
  }
}

typedef struct
{
  unsigned long remaining_length,
                remaining_symbols;
} MultisetItem;

DECLAREARRAYSTRUCT(MultisetItem);

#define MULTISET_PUSH_STACK(TYPE,ELEM)\
        STOREINARRAY(&stack,TYPE,256UL + 0.2 * stack.allocated##TYPE,ELEM)

uint8_t *multisets_list(unsigned long alphasize,unsigned long qgram_length)
{
  const unsigned long numelems = number_of_multisets(alphasize,qgram_length);
  uint8_t *multisets_space = gt_calloc(qgram_length * numelems,
                                       sizeof *multisets_space),
          *ms_next = multisets_space,
          ms_buf[MAX_QGRAM_LENGTH] = {0};
  ArrayMultisetItem stack;
  MultisetItem new_item;
#ifndef NDEBUG
  MultisetsInfo *ms_info = NULL;
  unsigned long a_idx = 0;

  if (qgram_length == 2 || qgram_length == 3)
  {
    ms_info = multisets_info_new(alphasize);
  }
#endif
  new_item.remaining_length = qgram_length;
  new_item.remaining_symbols = alphasize;
  INITARRAY(&stack,MultisetItem);
  MULTISET_PUSH_STACK(MultisetItem,new_item);
  while (stack.nextfreeMultisetItem > 0)
  {
    MultisetItem curr = stack.spaceMultisetItem[--stack.nextfreeMultisetItem];
    if (curr.remaining_length == 0)
    {
#ifndef NDEBUG
      if (ms_info != NULL)
      {
        unsigned long ms_encoded
          = multisets_encode(ms_info,ms_buf,qgram_length);
        assert(a_idx == ms_encoded);
        a_idx++;
      }
#endif
      memcpy(ms_next,ms_buf,sizeof *ms_next * qgram_length);
      ms_next += qgram_length;
    } else
    {
      if (curr.remaining_symbols > 0)
      {
        unsigned long idx, rest_of_alpha;
        uint8_t *ms_buf_ptr = ms_buf + qgram_length - curr.remaining_length;

        new_item.remaining_length = curr.remaining_length;
        new_item.remaining_symbols = curr.remaining_symbols - 1;
        MULTISET_PUSH_STACK(MultisetItem,new_item);
        for (idx = 1; idx <= curr.remaining_length; idx++)
        {
          new_item.remaining_length = curr.remaining_length - idx;
          new_item.remaining_symbols = curr.remaining_symbols - 1;
          MULTISET_PUSH_STACK(MultisetItem,new_item);
        }
        assert(alphasize >= curr.remaining_symbols);
        rest_of_alpha = alphasize - curr.remaining_symbols;
        for (idx = 0; idx < curr.remaining_length; idx++)
        {
          ms_buf_ptr[idx] = rest_of_alpha;
        }
      }
    }
  }
  assert (ms_next == multisets_space + qgram_length * numelems);
  FREEARRAY(&stack,MultisetItem);
#ifndef NDEBUG
  multisets_info_delete(ms_info);
#endif
  return multisets_space;
}
