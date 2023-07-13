#ifndef MULTISETS_H
#define MULTISETS_H
#include <stddef.h>
#include <inttypes.h>
#include <stdbool.h>
size_t number_of_multisets(size_t alphasize,size_t qgram_length);
uint8_t *multisets_list(unsigned long alphasize,unsigned long qgram_length);

typedef struct MultisetsInfo MultisetsInfo;

MultisetsInfo *multisets_info_new(unsigned long alphasize);
unsigned long multisets_encode(const MultisetsInfo *ms_info,const uint8_t *ms,
                               unsigned long qgram_length);
unsigned int **multisets_weights(const MultisetsInfo *ms_info,
                                 unsigned long qgram_length);
void multisets_info_delete(MultisetsInfo *ms_info);
#endif
