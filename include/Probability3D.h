/** \file IMP/domino3/Probability3D.h
    \brief A node for a single distance restraint.
 */

#ifndef IMPDOMINO3_PROBABILITY3D_H
#define IMPDOMINO3_PROBABILITY3D_H

IMPDOMINO3_BEGIN_NAMESPACE

void *memalign(size_t boundary, size_t size)
{
    void *pointer;
    int code = posix_memalign(&pointer,boundary,size);
    if (code != 0)
    {
        std::cerr<<"Error in memalign: Could not allocate memory by memalign. Please report this bug to developers\n";
        exit(3);
    }
    return pointer;
}

/** A node for a single distance restraint. */
class Probability3D {
private:
    int z_real_length,real_length;
    
    inline int getIndex(int x,int y,int z){
        IMP_USAGE_CHECK(x < x_length, "X is bigger than range");
        IMP_USAGE_CHECK(y < y_length, "Y is bigger than range");
        IMP_USAGE_CHECK(z < z_length, "Z is bigger than range");
        int index = x * y_length * z_real_length + y * z_real_length + z;
        return index;
    }
public:
    IMP::domino3::FP * array;
    int x_length,y_length,z_length,length;
    
    Probability3D(int x,int y, int z):
    x_length(x),
    y_length(y),
    z_length(z),
    length(x*y*z)
    {
        z_real_length=((z/4)+1)*4;
        real_length = (x_length*y_length*z_real_length);
        array = (IMP::domino3::FP *) memalign(16,(real_length) * sizeof(IMP::domino3::FP));
        std::fill (array,array+real_length,LogMathFunctions::convert_to_space(1));
    }
    
    ~Probability3D(){
        delete [] array;
    }
    
    inline IMP::domino3::FP get(int x,int y,int z){
        return  *(array+getIndex(x,y,z));
    }
    
    inline void set(int x,int y,int z,IMP::domino3::FP value){
        *(array + getIndex(x,y,z))  = value;
    }
    
    inline IMP::domino3::FP* get_ptr(int x,int y,int z){
        IMP::domino3::FP* ret_ptr =   array+( getIndex(x,y,z) );
        IMP_USAGE_CHECK(ret_ptr < array+real_length, "Address is out of ragen");

        return ret_ptr;
    }
    
    void show(){
        for(int x = 0; x < x_length; x++){
            std::cout << x << std::endl;
            for(int y = 0; y < y_length; y++){
                for(int z = 0; z < z_length; z++){
                    std::cout << get(x,y,z) << " " ;
                }
                std::cout << std::endl;
                
            }
        }
    }
    
    inline void normalize(){
        IMP::domino3::FP total = 0;
        for(int x = 0; x < x_length; x++){
            IMP::domino3::FP * it = get_ptr(x,0,0);
            total = std::accumulate(it, it + (y_length*z_real_length), 0.0);
            for(int y = 0; y < y_length; y++){
                for(int z = 0; z < z_length; z++){
                    set(x,y,z,IMP::domino3::LogMathFunctions::convert_to_space(get(x,y,z)/total));
                }
            }
        }
    }
    
};


IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_PROBABILITY3D_H */
