#ifndef IMPDOMINO3_MATH_H
#define IMPDOMINO3_MATH_H

IMPDOMINO3_BEGIN_NAMESPACE

class MathFunctions {
public:
	virtual double add(double val1, double val2)=0;

	virtual double mult(double val1,double val2)=0;

	virtual double convert_to_space(double val1)=0;

        virtual void normalize(double * it, unsigned int size)=0;

        virtual double sum(double * vals, unsigned int size)=0;

        virtual double convert_to_linear(double val1)=0;
};

IMPDOMINO3_END_NAMESPACE
#endif // IMPDOMINO3_NODE_H
