#include <IMP/domino3/LinearMathFunctions.h>
IMPDOMINO3_BEGIN_NAMESPACE
inline double LinearMathFunctions::add(double val1, double val2){
        return val1+val2;
}

inline double LinearMathFunctions::mult(double val1,double val2){
        return val1*val2;
}

inline double LinearMathFunctions::dev(double val1,double val2){
        return val1/val2;
}

inline double LinearMathFunctions::convert_to_space(double val1){
        return val1;
}

void LinearMathFunctions::normalize(double * it,
                                    unsigned int size) {
        double total = std::accumulate(it, it + size, 0.0);
//        IMP_USAGE_CHECK(total > .001,
//                        "Total is too small to be reliable: " << total);
        for (unsigned int i = 0; i < size; ++i) {
            it[i] /= total;
        }
}

inline double LinearMathFunctions::convert_to_linear(double val1){
	return val1;
}

inline double LinearMathFunctions::sum(double * vals, unsigned int size) {
	return std::accumulate(vals, vals + size, 0.0);
}

IMPDOMINO3_END_NAMESPACE
