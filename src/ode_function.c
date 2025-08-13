#include <math.h>
#include "ode_function.h"

double ode_function(double x, double y) {
    return 10 * (12 * cos(2 * M_PI * x) - y);
}