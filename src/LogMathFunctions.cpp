#include <IMP/domino3/LogMathFunctions.h>

IMPDOMINO3_BEGIN_NAMESPACE
// Impl. of log exp trick
// http://machineintelligence.tumblr.com/post/4998477107/the-log-sum-exp-trick
double LogMathFunctions::add(double val1, double val2){
        double max=std::max(val1,val2);
        double sum = 0;
        sum += exp(val1 - max);
        sum += exp(val2 - max);
        return max+log(sum);
}

double LogMathFunctions::dev(double val1,double val2){
        return val1-val2;
}

double LogMathFunctions::sum(double * vals, unsigned int size){
    double max=*std::max_element(vals,vals+size);
    double sum = 0;
    for(int i = 0; i < size; i++){
        sum += exp(vals[i] - max);
    }
    return max+log(sum);
}


double LogMathFunctions::mult(double val1,double val2){
        return val1+val2;
}

double LogMathFunctions::convert_to_space(double val1){
        return (val1==0) ? -DBL_MAX : log(val1);
}

double LogMathFunctions::convert_to_linear(double val1){
	return exp(val1);
}

void LogMathFunctions::normalize(double * it,
                                    unsigned int size) {
    double sum_vec=sum(it,size);
    for(int i = 0; i < size; i++){
        it[i]-=sum_vec;
    }
}
IMPDOMINO3_END_NAMESPACE
