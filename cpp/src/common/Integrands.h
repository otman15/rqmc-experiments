#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include <vector>

double sumueu(const std::vector<double>& x);

double smoothperb4(
    const std::vector<double>& x,
    double expw);

double mc2(const std::vector<double>& x);

double polynomial(const std::vector<double>& x);

double oscillatory(const std::vector<double>& x);

double gaussian(const std::vector<double>& x);

double smoothgauss(const std::vector<double>& x);

double piecelingauss(const std::vector<double>& x);

double indsumnormal(const std::vector<double>& x);

#endif
