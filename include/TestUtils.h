#pragma once

#include "VRandGen.h"

inline const char* queryModeName(VRandGenQueryMode qm)
{
    switch (qm) {
        case QM_Any: return "AnySize";
        case QM_Scalar: return "Scalar";
        case QM_Block16: return "Block16";
        case QM_StateSize: return "State";
        default: THROW("how did we get here?")
    }
}
