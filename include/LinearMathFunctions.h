#ifndef IMPDOMINO3_LINEARMATH_H
#define IMPDOMINO3_LINEARMATH_H 
#include <IMP/domino3/MathFunctions.h>
IMPDOMINO3_BEGIN_NAMESPACE

/** A node for a single distance restarint. */
class IMPDOMINO3EXPORT LinearMathFunctions: public MathFunctions {
 public:
        virtual double add(double val1, double val2) IMP_OVERRIDE;

        virtual double mult(double val1,double val2) IMP_OVERRIDE;

        virtual double dev(double val1,double val2) IMP_OVERRIDE;
    
        virtual double convert_to_space(double val1) IMP_OVERRIDE;

        virtual void normalize(double * it, unsigned int size) IMP_OVERRIDE;

        virtual double sum(double * vals, unsigned int size) IMP_OVERRIDE;

        virtual double convert_to_linear(double val1) IMP_OVERRIDE;
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_DISTANCE_NODE_H
