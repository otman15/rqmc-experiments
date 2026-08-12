#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include <vector>

double sumueu(const std::vector<double>& x);

double smoothperb4(
    const std::vector<double>& x,
    double expw);

double mc2(const std::vector<double>& x);

double polynomial(const std::vector<double>& x);

#endif
