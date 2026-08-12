#include <iostream>

#include "BurleyQuantLib.h"

int main()
{
    unsigned int s = 2;
    unsigned int k = 3;
    long seed = 12345;

    BurleyQuantLib sequence(s, k, seed);

    for (std::uint64_t i = 0; i < sequence.numPoints(); ++i)
    {
        for (unsigned int j = 0; j < s; ++j)
        {
            std::cout << sequence.coordinate(j) << " ";
        }

        std::cout << '\n';
        sequence.next();
    }

    return 0;
}
