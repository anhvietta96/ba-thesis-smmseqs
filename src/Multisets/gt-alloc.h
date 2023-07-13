#ifndef GT_ALLOC_H
#define GT_ALLOC_H
#include <stdlib.h>

static inline void *gt_malloc_generic(
       __attribute__ ((unused)) const char *filename,
       __attribute__ ((unused)) int line,
       size_t size)
{
  /*printf("%s,%s,%d\n",__func__,filename,line);*/
  return malloc(size);
}

static inline void *gt_calloc_generic(
        __attribute__ ((unused)) const char *filename,
        __attribute__ ((unused)) int line,
        size_t nmemb,size_t size)
{
  /*printf("%s,%s,%d\n",__func__,filename,line);*/
  return calloc(nmemb,size);
}

static inline void *gt_realloc_generic(
         __attribute__ ((unused)) const char *filename,
         __attribute__ ((unused)) int line,
        void *ptr,size_t size)
{
  /*printf("%s,%s,%d\n",__func__,filename,line);*/
  return realloc(ptr,size);
}

#define gt_malloc(X)    gt_malloc_generic(__FILE__,__LINE__,X)
#define gt_realloc(X,Y) gt_realloc_generic(__FILE__,__LINE__,X,Y)
#define gt_calloc(X,Y)  gt_calloc_generic(__FILE__,__LINE__,X,Y)
#define gt_free(X)      free(X)

#endif
