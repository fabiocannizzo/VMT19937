#pragma once

#include "RandGen.h"
#include "cpu.h"

enum QryMode { QM_Scalar, QM_Block16, QM_Any };

inline const char* queryModeName(QryMode qm)
{
    switch (qm) {
        case QM_Any: return "AnySize";
        case QM_Scalar: return "Scalar";
        case QM_Block16: return "Block16";
        default: THROW("how did we get here?")
    }
}

