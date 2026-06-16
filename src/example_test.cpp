/*
Posc++ - Library for Power-Series Composition
The files in this project are maintained in the GitHub repository Posc++, available at
https://github.com/V3du4rd0/posc-.
*/

#include "PowerSeries.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <map>
#include <sys/stat.h>
#include <filesystem>
#include <cstdio>

#include "functionENV.h"



int main(int argc, char *argv[]) {
  std::cout << "Test! " << std::endl;

// Sin( 2 Cos( X ) )
cstm_float_t u0 =  val(0.25)*Pi; // evaluation point fot Taylor series
const int N = 6; // number of terms
auto u = FE::Variable(u0, N); // independent variable

// custom function sin2cos(u)	
auto sin2cos = FE::SIN_K(FE::CTE_MULT(val(2.0), FE::COS_K(u, N)  ,N),N);
// retrieve coefficients
std::vector<cstm_float_t> coeffs_STA = sin2cos->get_Taylor_coefficients(N);

std::cout<<"sin_2cos at u0 = "; print_result(u0); std::cout<<"\n\n";
for(int i=0; i<N; i++){
        std::cout<< i<<"-th term "; print_result(coeffs_STA[i]); std::cout<<"\n";
}

return 0;
}
