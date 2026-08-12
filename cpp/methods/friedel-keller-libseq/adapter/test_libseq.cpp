#include <iostream>
#include <cstdlib>

#include "RandomizedTSSequence.h"

int main()
{
    UL_int dimension = 2;
    UL_int length = 8;

    srand48(12345);

    Ring R(2);
    C_matrix C(&R, "Niederreiter", dimension);

    RandomizedTSSequence_base_2 seq(&C, dimension, length);

    for (UL_int i = 0; i < length; ++i)
    {
        for (UL_int j = 0; j < dimension; ++j)
        {
            std::cout << seq[j] << " ";
        }

        std::cout << std::endl;

        if (i + 1 < length)
            ++seq;
    }

    return 0;
}
