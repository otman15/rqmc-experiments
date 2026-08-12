#include "RqmcReplications.h"

#include "Integrands.h"

#include <stdexcept>

double evaluateRqmcModel(
    const std::string& model,
    const std::vector<double>& x)
{
    if (model == "sumueu")
        return sumueu(x);

    if (model == "smoothperb4")
        return smoothperb4(x, 1.0);

    if (model == "mc2")
        return mc2(x);

    if (model == "polynomial")
        return polynomial(x);

    throw std::invalid_argument(
        "Unknown RQMC model: " + model);
}
