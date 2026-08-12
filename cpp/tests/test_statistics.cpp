#include <iostream>
#include <vector>

#include "Statistics.h"

int main()
{
    const std::vector<double> odd = {
        5.0, 1.0, 9.0, 3.0, 7.0
    };

    const std::vector<double> even = {
        8.0, 2.0, 6.0, 4.0
    };

    std::cout << "odd median  = " << median(odd) << '\n';
    std::cout << "even median = " << median(even) << '\n';

    return 0;
}
