#include "nui.h"

std::ostream& operator<<(std::ostream& os, const nui& per) {
    os << per.a;
    return os;
}