#ifndef RAMSEY_FIXED_P_CUH_INCLUDED
#define RAMSEY_FIXED_P_CUH_INLCUDED

__device__ void calculate_ramsey_fixed_p( 
	float const * parameters,
    int const n_fits,
    int const n_points,
    float * value,
    float * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
	float * user_info_float = (float*)user_info;
    float x = 0.0f;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(float) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(float) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }
	
	//////////////////////////// values //////////////////////////////
	
	// parameters: [A1 A2 c f1 f2 t2star x1 x2] exp(-(point_index./t2star)^1)*(A1*cos(2*pi*f1*(point_index - x1)) + A2*cos(2*pi*f2*(point_index-x2))) + c
	
	float const * p = parameters;
	
	float const pi = 3.14159;
	float const t2arg = x/p[5];
	float const ex = exp(-t2arg);
	float const phasearg1 = 2*pi*p[3]*(x - p[6]);
	float const phasearg2 = 2*pi*p[4]*(x - p[7]);
	float const cos1 =  cos(phasearg1);
	float const sin1 =  sin(phasearg1);
	float const cos2 =  cos(phasearg2);
	float const sin2 =  sin(phasearg2);
	
	value[point_index] = ex*(p[0]*cos1 + p[1]*cos2) + p[2] ; // formula calculating fit model values
	
	/////////////////////////// derivatives ///////////////////////////
	float * current_derivative = derivative + point_index;
	
	current_derivative[0 * n_points ] = ex*cos1 ; // formula calculating derivative values with respect to parameters[0]
	current_derivative[1 * n_points ] = ex*cos2;
	current_derivative[2 * n_points ] = 1.f;
	current_derivative[3 * n_points ] = -p[0]*2*pi*(x-p[6])*ex*sin1;
	current_derivative[4 * n_points ] = -p[1]*2*pi*(x-p[7])*ex*sin2;
	current_derivative[5 * n_points ] = 1.f/(p[5]*p[5])*x*ex*(p[0]*cos1+p[1]*cos2);
	current_derivative[6 * n_points ] = p[0]*2*pi*p[3]*sin1*ex;
	current_derivative[7 * n_points ] = p[1]*2*pi*p[4]*sin2*ex;
	
	
}

#endif