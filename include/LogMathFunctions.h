/** \file IMP/domino3/LogMathFunctions.h
 */

#ifndef IMPDOMINO3_LOGMATHFUNCTIONS_H
#define IMPDOMINO3_LOGMATHFUNCTIONS_H

#include <xmmintrin.h>

IMPDOMINO3_BEGIN_NAMESPACE
typedef float FP;
/** A node for a single distance restraint. */
class LogMathFunctions {
public:
    
    static inline float get_first (const __m128 vec){return _mm_cvtss_f32(_mm_shuffle_ps(vec,vec, _MM_SHUFFLE2(0,0)));}
    static inline float get_second(const __m128 vec){return _mm_cvtss_f32(_mm_shuffle_ps(vec,vec, _MM_SHUFFLE2(0,1)));}
    static inline float get_third (const __m128 vec){return _mm_cvtss_f32(_mm_shuffle_ps(vec,vec, _MM_SHUFFLE2(1,0)));}
    static inline float get_fourth(const __m128 vec){return _mm_cvtss_f32(_mm_shuffle_ps(vec,vec, _MM_SHUFFLE2(1,1)));}
    // Impl. of log exp trick
    // http://machineintelligence.tumblr.com/post/4998477107/the-log-sum-exp-trick
    static inline FP add(FP val1, FP val2){
        FP max=std::max(val1,val2);
        max=std::max(1.0f,max);
        
        FP sum = 0;
        
        sum += convert_to_linear(val1 - max);
        sum += convert_to_linear(val2 - max);
        FP ret_val = max+convert_to_space(sum);
        return ret_val;
    }
    
    
    static inline __m128 mult_sse(__m128 term1,__m128 term2){
        return _mm_add_ps(term1,term2);
    }
    
    
    static inline __m128 add_sse(__m128 val1, __m128 val2){
        __m128 one =  _mm_set1_ps(1.0f);
        __m128 max = _mm_max_ps(val1,val2);
        max = _mm_max_ps(one,max);
        __m128 sum =  _mm_set1_ps(0.0f);
        __m128 val1_minus_max=_mm_fpow2(_mm_sub_ps(val1, max));
        sum = _mm_add_ps(sum,val1_minus_max);
        __m128 val2_minus_max=_mm_fpow2(_mm_sub_ps(val2, max));
        sum = _mm_add_ps(sum,val2_minus_max);
        __m128 ret_val=_mm_add_ps(max,_mm_flog2_ps(sum));
        return ret_val;
    }
    
    
    static inline FP dev(FP val1,FP val2){
        return val1-val2;
    }
    
    static inline FP sum(FP * vals, unsigned int size){
        FP max=*std::max_element(vals,vals+size);
        max=std::max(1.0f,max);
        FP sum = 0;
        for(unsigned i = 0; i < size; i++){
            sum += convert_to_linear(vals[i] - max);
        }
        return max+convert_to_space(sum);
    }
    
    
    static inline FP mult(FP val1,FP val2){
        return val1+val2;
    }
    
    static inline FP convert_to_space(FP val1){
        return flog2(val1);
    }
    
    static inline FP convert_to_linear(FP val1){
        return powf(2, val1);
    }
    
    static void normalize(FP * it,
                          unsigned int size) {
        FP sum_vec=sum(it,size);
        for(unsigned i = 0; i < size; i++){
            it[i]-=sum_vec;
        }
    }
    

    
private:
    /////////////////////////////////////////////////////////////////////////////////////
    // fast log base 2
    /////////////////////////////////////////////////////////////////////////////////////
    // Fast log2
    // ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
    // Maximum deviation: +/- 2.1E-5
    // Run time: ~1.2E-8s on Intel core2 2.13GHz, log2(): 5.4E-8s
    // For a negative argument, -128 is returned.
    // The function makes use of the representation of 4-byte floating point numbers:
    // seee eeee emmm mmmm mmmm mmmm mmmm mmmm
    // s is the sign, eee eee e gives the exponent + 127 (in hex: 0x7f).
    // The following 23 bits give the mantisse, the binary digits after the decimal
    // point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
    // Therefore,  log2(x) = eeeeeeee-127 + log2(1.mmmmmm...)
    //                     = eeeeeeee-127 + log2(1+y),  where y = 0.mmmmmm...
    //                     ~ eeeeeeee-127 + ((a*y+b)*y+c)*y
    // The coefficients a, b  were determined by a least squares fit, and c=1-a-b to get 1 at y=1.
    // Lower/higher order polynomials may be used for faster or more precise calculation:
    // Order 1: log2(1+y) ~ y
    // Order 2: log2(1+y) = (a*y + 1-a)*y, a=-0.3427
    //  => max dev = +/- 8E-3, run time ~ ?
    // Order 3: log2(1+y) = ((a*y+b)*y + 1-a-b)*y, a=0.1564, b=-0.5773
    //  => max dev = +/- 1E-3, run time ~ ?
    // Order 4: log2(1+y) = (((a*y+b)*y+c)*y + 1-a-b-c)*y, a=-0.0803 b=0.3170 c=-0.6748
    //  => max dev = +/- 1.4E-4, run time ~ ?
    // Order 5: log2(1+y) = ((((a*y+b)*y+c)*y+d)*y + 1-a-b-c-d)*y,
    //     a=0.0440047 b=-0.1903190 c=0.4123442 d=-0.7077702 1-a-b-c-d=1.441740
    //  => max dev = +/- 2.1E-5, run time ~ 1.2E-8s
    static inline float flog2(float x) {
        if (x <= 0)
            return -128;
        int *px = (int*) (&x);      // store address of float as pointer to long int
        float e = (float) (((*px & 0x7F800000) >> 23) - 0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
        *px = ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0)
        x -= 1.0;         // and calculate x-1.0
        x *=
        (1.441740
         + x
         * (-0.7077702
            + x
            * (0.4123442
               + x
               * (-0.1903190
                  + x
                  * 0.0440047)))); // 5'th order polynomial approx. of log(1+x)
        return x + e;
    }
    
    // Fast SSE2 log2 for four floats
    // Calculate integer of log2 for four floats in parallel with SSE2
    // Maximum deviation: +/- 2.1E-5
    // Run time: ~5.6ns on Intel core2 2.13GHz.
    // For a negative argument, nonsense is returned. Otherwise, when <1E-38, a value
    // close to -126 is returned and when >1.7E38, +128 is returned.
    // The function makes use of the representation of 4-byte floating point numbers:
    // seee eeee emmm mmmm mmmm mmmm mmmm mmmm
    // s is the sign, eee eee e gives the exponent + 127 (in hex: 0x7f).
    // The following 23 bits give the mantisse, the binary digits after the decimal
    // point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
    // Therefore,  log2(x) = eeeeeeee-127 + log2(1.mmmmmm...)
    //                     = eeeeeeee-127 + log2(1+y),  where y = 0.mmmmmm...
    //                     ~ eeeeeeee-127 + ((a*y+b)*y+c)*y
    // The coefficients a, b  were determined by a least squares fit, and c=1-a-b to get 1 at y=1.
    // Lower/higher order polynomials may be used for faster or more precise calculation:
    // Order 1: log2(1+y) ~ y
    // Order 2: log2(1+y) = (a*y + 1-a)*y, a=-0.3427
    //  => max dev = +/- 8E-3, run time ~ 3.8ns
    // Order 3: log2(1+y) = ((a*y+b)*y + 1-a-b)*y, a=0.1564, b=-0.5773
    //  => max dev = +/- 1E-3, run time ~ 4.4ns
    // Order 4: log2(1+y) = (((a*y+b)*y+c)*y + 1-a-b-c)*y, a=-0.0803 b=0.3170 c=-0.6748
    //  => max dev = +/- 1.4E-4, run time ~ 5.0ns?
    // Order 5: log2(1+y) = ((((a*y+b)*y+c)*y+d)*y + 1-a-b-c-d)*y, a=0.0440047 b=-0.1903190 c=0.4123442 d=-0.7077702
    //  => max dev = +/- 2.1E-5, run time ~ 5.6ns?
    
    //#ifdef HH_SSE2
    static inline __m128 _mm_flog2_ps(__m128 X)
    //inline __m128 logariSSE(__m128 X)
    {
        const __m128i CONST32_0x7f = _mm_set_epi32(0x7f, 0x7f, 0x7f, 0x7f);
        const __m128i CONST32_0x7fffff = _mm_set_epi32(0x7fffff, 0x7fffff, 0x7fffff,
                                                       0x7fffff);
        const __m128i CONST32_0x3f800000 = _mm_set_epi32(0x3f800000, 0x3f800000,
                                                         0x3f800000, 0x3f800000);
        const __m128 CONST32_1f = _mm_set_ps(1.0, 1.0, 1.0, 1.0);
        // const float a=0.1564, b=-0.5773, c=1.0-a-b;  // third order
        const float a = 0.0440047, b = -0.1903190, c = 0.4123442, d = -0.7077702,
        e = 1.0 - a - b - c - d; // fifth order
        const __m128 CONST32_A = _mm_set_ps(a, a, a, a);
        const __m128 CONST32_B = _mm_set_ps(b, b, b, b);
        const __m128 CONST32_C = _mm_set_ps(c, c, c, c);
        const __m128 CONST32_D = _mm_set_ps(d, d, d, d);
        const __m128 CONST32_E = _mm_set_ps(e, e, e, e);
        __m128i E; // exponents of X
        __m128 R; //  result
        
        E = _mm_srli_epi32((__m128i) X, 23); // shift right by 23 bits to obtain exponent+127
        E = _mm_sub_epi32(E, CONST32_0x7f); // subtract 127 = 0x7f
        X = (__m128) _mm_and_si128((__m128i) X, CONST32_0x7fffff); // mask out exponent => mantisse
        X = (__m128) _mm_or_si128((__m128i) X, CONST32_0x3f800000); // set exponent to 127 (i.e., 0)
        X = _mm_sub_ps(X, CONST32_1f); // subtract one from mantisse
        R = _mm_mul_ps(X, CONST32_A); // R = a*X
        R = _mm_add_ps(R, CONST32_B); // R = a*X+b
        R = _mm_mul_ps(R, X); // R = (a*X+b)*X
        R = _mm_add_ps(R, CONST32_C); // R = (a*X+b)*X+c
        R = _mm_mul_ps(R, X); // R = ((a*X+b)*X+c)*X
        R = _mm_add_ps(R, CONST32_D); // R = ((a*X+b)*X+c)*X+d
        R = _mm_mul_ps(R, X); // R = (((a*X+b)*X+c)*X+d)*X
        R = _mm_add_ps(R, CONST32_E); // R = (((a*X+b)*X+c)*X+d)*X+e
        R = _mm_mul_ps(R, X); // R = ((((a*X+b)*X+c)*X+d)*X+e)*X ~ log2(1+X) !!
        R = _mm_add_ps(R, _mm_cvtepi32_ps(E)); // convert integer exponent to float and add to mantisse
        return R;
    }
    //#endif
    
    // This function returns log2 with a max absolute deviation of +/- 1.5E-5 (typically 0.8E-5).
    // It takes 0.80E-8 s  whereas log2(x) takes 5.4E-7 s. It is hence 9.4 times faster.
    // It makes use of the representation of 4-byte floating point numbers:
    // seee eeee emmm mmmm mmmm mmmm mmmm mmmm
    // s is the sign,
    // the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
    // The following 23 bits give the mantisse, the binary digits after the decimal
    // point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
    // In the code, *(int *)&x is an integer which contains the bytes as the
    // floating point variable x is represented in memory. The expression
    //     (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee,
    // i.e., the largest integer that is smaller than log2(x) (e.g. -1 for 0.9).
    static inline float fast_log2(float x) {
        
        union { float f; uint32_t i; } vx = { x };
        union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
        float y = vx.i;
        y *= 1.1920928955078125e-7f;
        
        return (y - 124.22551499f
                - 1.498030302f * mx.f
                - 1.72587999f / (0.3520887068f + mx.f));
    }
    
    

    
    /////////////////////////////////////////////////////////////////////////////////////
    // fast 2^x
    // ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
    // Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
    // Speed: 2.1E-8s (2.3E-8s) per call! (exp(): 8.5E-8, pow(): 1.7E-7)
    // Internal representation of float number according to IEEE 754:
    //   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
    //                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000
    //   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
    /////////////////////////////////////////////////////////////////////////////////////
    static inline float fpow2(float x) {
        if (x > FLT_MAX_EXP)
            return FLT_MAX;
        if (x < FLT_MIN_EXP)
            return 0.0f;
        int *px = (int*) (&x);      // store address of float as pointer to long int
        float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        int lx = *((int*) &tx) - 0x4b400000;   // integer value of x
        float dx = x - (float) (lx);             // float remainder of x
        //   x = 1.0f + dx*(0.69606564f           // cubic apporoximation of 2^x for x in the range [0, 1]
        //            + dx*(0.22449433f           // Gives relative deviation < 1.5E-4
        //            + dx*(0.07944023f)));       // Speed: 1.9E-8s
        x = 1.0f + dx * (0.693019f // polynomial approximation of 2^x for x in the range [0, 1]
                         + dx * (0.241404f             // Gives relative deviation < 4.6E-6
                                 + dx * (0.0520749f            // Speed: 2.1E-8s
                                         + dx * 0.0134929f)));
        //   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
        //            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
        //            + dx*(0.0558282f            // Speed: 2.3E-8s
        //            + dx*(0.00898898f
        //            + dx* 0.00187682f ))));
        *px += (lx << 23);                     // add integer power of 2 to exponent
        return x;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    // SSE 2^x for four floats
    // Calculate float of 2pow(x) for four floats in parallel with SSE2
    // ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
    // Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
    //
    // Internal representation of float number according to IEEE 754 (__m128 --> 4x):
    //   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
    //                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000
    //   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
    /////////////////////////////////////////////////////////////////////////////////////
    //#ifdef HH_SSE2
    static inline __m128 _mm_fpow2(__m128 X) {
        
        __m128i* xPtr = (__m128i*) &X;	// store address of float as pointer to int
        
        const __m128 CONST32_05f = _mm_set1_ps(0.5f); // Initialize a vector (4x32) with 0.5f
        // (3 << 22) --> Initialize a large integer vector (shift left)
        const __m128i CONST32_3i = _mm_set_epi32(3, 3, 3, 3);
        const __m128i CONST32_3shift22 = _mm_slli_epi32(CONST32_3i, 22);
        const __m128 CONST32_1f = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f);
        const __m128 CONST32_FLTMAXEXP = _mm_set1_ps(FLT_MAX_EXP);
        const __m128 CONST32_FLTMAX = _mm_set1_ps(FLT_MAX);
        const __m128 CONST32_FLTMINEXP = _mm_set1_ps(FLT_MIN_EXP);
        // fifth order
        const __m128 CONST32_A = _mm_set1_ps(0.00187682f);
        const __m128 CONST32_B = _mm_set1_ps(0.00898898f);
        const __m128 CONST32_C = _mm_set1_ps(0.0558282f);
        const __m128 CONST32_D = _mm_set1_ps(0.240153f);
        const __m128 CONST32_E = _mm_set1_ps(0.693153f);
        
        __m128 tx;
        __m128i lx;
        __m128 dx;
        __m128 result = _mm_set1_ps(0.0f);
        __m128 maskedMax = _mm_set1_ps(0.0f);
        __m128 maskedMin = _mm_set1_ps(0.0f);
        
        // Check wheter one of the values is bigger or smaller than FLT_MIN_EXP or FLT_MAX_EXP
        // The correct FLT_MAX_EXP value is written to the right place
        maskedMax = _mm_cmpgt_ps(X, CONST32_FLTMAXEXP);
        maskedMin = _mm_cmpgt_ps(X, CONST32_FLTMINEXP);
        maskedMin = _mm_xor_ps(maskedMin, maskedMax);
        // If a value is bigger than FLT_MAX_EXP --> replace the later result with FLTMAX
        maskedMax = _mm_and_ps(CONST32_FLTMAX, _mm_cmpgt_ps(X, CONST32_FLTMAXEXP));
        
        tx = _mm_add_ps((__m128 ) CONST32_3shift22, _mm_sub_ps(X, CONST32_05f)); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        
        lx = _mm_cvtps_epi32(tx);										// integer value of x
        
        dx = _mm_sub_ps(X, _mm_cvtepi32_ps(lx));						// float remainder of x
        
        //   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
        //            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
        //            + dx*(0.0558282f            // Speed: 2.3E-8s
        //            + dx*(0.00898898f
        //            + dx* 0.00187682f ))));
        X = _mm_mul_ps(dx, CONST32_A);
        X = _mm_add_ps(CONST32_B, X);	// add constant B
        X = _mm_mul_ps(dx, X);
        X = _mm_add_ps(CONST32_C, X);	// add constant C
        X = _mm_mul_ps(dx, X);
        X = _mm_add_ps(CONST32_D, X);	// add constant D
        X = _mm_mul_ps(dx, X);
        X = _mm_add_ps(CONST32_E, X);	// add constant E
        X = _mm_mul_ps(dx, X);
        X = _mm_add_ps(X, CONST32_1f);	// add 1.0f
        
        __m128i lxExp = _mm_slli_epi32(lx, 23); // add integer power of 2 to exponent
        
        *xPtr = _mm_add_epi32(*xPtr, lxExp); // add integer power of 2 to exponent
        
        // Add all Values that are greater than min and less than max
        result = _mm_and_ps(maskedMin, X);
        // Add MAX_FLT values where entry values were > FLT_MAX_EXP
        result = _mm_or_ps(result, maskedMax);
        
        return result;
    }
    
    //#endif
};



IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_LOGMATHFUNCTIONS_H */
