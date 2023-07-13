#ifndef MINMAX_H
#define MINMAX_H
#ifndef MAX
#define MAX(X,Y) ((X) < (Y) ? (Y) : (X))
#endif
#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX3
#define MAX3(X, Y, Z) (((X)>(Y))?((X)>(Z)?(X):(Z)):((Y)>(Z)?(Y):(Z)))
#endif
#ifndef MIN3
#define MIN3(X, Y, Z) (((X)<(Y))?((X)<(Z)?(X):(Z)):((Y)<(Z)?(Y):(Z)))
#endif
#endif
