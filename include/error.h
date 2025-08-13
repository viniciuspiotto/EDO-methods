#ifndef ERROR_H
#define ERROR_H
#include "methods.h"

double * error (Result res, Result correct);
double iterationError (Result res, Result correct, int start, int end);

#endif