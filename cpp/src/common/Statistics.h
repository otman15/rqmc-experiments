#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>

double mean(const std::vector<double>& values);

double median(const std::vector<double>& values);

double sampleVariance(const std::vector<double>& values);

double skewness(const std::vector<double>& values);

double excessKurtosis(const std::vector<double>& values);

#endif
