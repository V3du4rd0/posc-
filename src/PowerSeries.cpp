/*
Posc++ - Library for Power-Series Composition
The files in this project are maintained in the GitHub repository Posc++, available at
https://github.com/V3du4rd0/posc-.
*/

#include "PowerSeries.h"
#include "functionENV.h"

template <typename T>
void print_result(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        char buffer[128];
        quadmath_snprintf(buffer, sizeof(buffer), "%.36Qg", value); // Adjust precision as needed
        std::cout << buffer ;//<< std::endl;
    } else {
        if constexpr (std::is_same<T, double>::value) {
            std::cout << std::fixed << std::setprecision(15) <<  value;// << std::endl;
        } else {
            std::cout << std::fixed << std::setprecision(7) <<  value;// << std::endl;
        }
    }
}

template void print_result(float);
template void print_result(double);
template void print_result(__float128);



cstm_float_t val(double value) {
    std::string str_value = std::to_string(value);

    if constexpr (std::is_same<cstm_float_t, __float128>::value) {
        return strtoflt128(str_value.c_str(), NULL);
    } else if constexpr (std::is_same<cstm_float_t, double>::value) {
        return std::strtod(str_value.c_str(), NULL);
    } else {
        return std::strtof(str_value.c_str(), NULL);
    }
}

// #

Function::Function(int N) : max_order(N) {}

Function::~Function() = default;

std::vector<cstm_float_t> Function::get_Taylor_coefficients(int N) const {
    std::vector<cstm_float_t> coeffs;
    for (int i = 0; i < N; ++i) {
        coeffs.push_back(get_coefficient(i));
    }
    return coeffs;
}


Variable::Variable(cstm_float_t p_, int N) : Function(N), p(p_) {
    coeffs.resize(N, 0.0);
    coeffs[0] = p;
    /*
    if (N >= 1){
        coeffs[1] = 1.0;
    }*/
    coeffs[1] = val(1.0);
    for(int i=2;i<N;i++){
        coeffs[i] = val(0.0);
    }
}

cstm_float_t Variable::get_coefficient(int k) const {
    //assert(k >= 0 && k < coeffs.size());
    return coeffs[k];
}


CTE::CTE(cstm_float_t val, int N) : Function(N), value(val) {
    coeffs.resize(N, 0.0);
    coeffs[0] = val;
}

cstm_float_t CTE::get_coefficient(int k) const {
        return (k == 0 ? value : 0.0);
}


POWER2::POWER2(std::shared_ptr<Function> f_, int N) : Function(N), f(f_) {
    cstm_float_t aux = f->get_coefficient(0);
    std::vector<cstm_float_t> pol;

    pol.push_back(aux * aux);
    for (int k = 1; k < N ; ++k) {
        cstm_float_t suma = val(0.0);
        for (int j = 0; j <= k; ++j) {
           suma += (f->get_coefficient(j))*(f->get_coefficient(k-j));
        }
        pol.push_back(suma);
    }
    this->coeffs = pol;
}

cstm_float_t POWER2::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);
        if (k == 0){
            return c_pow(f->get_coefficient(0),val(2.0));
        }
        return coeffs[k];
}

SIN_K::SIN_K(std::shared_ptr<Function> f_, int N) : Function(N), f(f_) {
    //cstm_float_t aux = f->get_coefficient(0);
    std::vector<cstm_float_t> f_Taylor = f->get_Taylor_coefficients(N);
    cstm_float_t aux = f_Taylor[0];

    //if (f_Taylor.size() < static_cast<size_t>(N + 1)) {
    //    f_Taylor.resize(N + 1, val(0.0));
    //}

        //C = new cstm_float_t[f->degree + 1];
        //S = new cstm_float_t[f->degree + 1];
    std::vector<cstm_float_t> C;
    std::vector<cstm_float_t> S;

    C.push_back(c_cos(aux));
        S.push_back(c_sin(aux));

        for (int k = 1; k < N ; ++k) {
            cstm_float_t termj_c = val(0.0);
            cstm_float_t termj_s = val(0.0);
            for (int j = 0; j < k; ++j) {
            //aux = (f->get_coefficient(k-j));//f->TAYLOR_K(k - j);
            aux = f_Taylor[k - j];

                termj_s += val(k - j) * C[j] * aux;
                termj_c += val(k - j) * S[j] * aux;
            }
            termj_c = val(-1.0)*termj_c / val(k);
            termj_s = termj_s / val(k);
            C.push_back(termj_c);
            S.push_back( termj_s);
        }

        this->coeffs = S;
}

cstm_float_t SIN_K::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);
        //if (k == 0){
        //    return c_sin(f->get_coefficient(0));
        //}
        //return coeffs[k];
        assert(k >= 0 && k < coeffs.size());
        return coeffs[k];
}



COS_K::COS_K(std::shared_ptr<Function> f_, int N) : Function(N), f(f_) {
    //cstm_float_t aux = f->get_coefficient(0);
    std::vector<cstm_float_t> f_Taylor = f->get_Taylor_coefficients(N);
    cstm_float_t aux = f_Taylor[0];

    std::vector<cstm_float_t> C;
    std::vector<cstm_float_t> S;

    C.push_back(c_cos(aux));
        S.push_back(c_sin(aux));

        for (int k = 1; k < N ; ++k) {
            cstm_float_t termj_c = val(0.0);
            cstm_float_t termj_s = val(0.0);
            for (int j = 0; j < k; ++j) {
            aux =f->get_coefficient(k - j);
                termj_s += val(k - j) * C[j] * aux;
                termj_c += val(k - j) * S[j] * aux;
            }
            termj_c = val(-1.0)*termj_c / val(k);
            termj_s = termj_s / val(k);
            C.push_back(termj_c);
            S.push_back( termj_s);
        }

        this->coeffs = C;

}

cstm_float_t COS_K::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);
        //if (k == 0){
        //    return c_cos(f->get_coefficient(0));
        //}
        //return coeffs[k];
        assert(k >= 0 && k < coeffs.size());
        return coeffs[k];
}


TAN_K::TAN_K(std::shared_ptr<Function> f_, int N) : Function(N), f(f_) {
    //cstm_float_t aux = f->get_coefficient(0);
    std::vector<cstm_float_t> f_coeffs = f_->get_Taylor_coefficients(N);
    cstm_float_t aux = f_coeffs[0];

    auto cos_ptr = std::make_shared<COS_K>(f_, N);
    auto c2 = FE::POWER2(cos_ptr, N);
    //auto c2 = FE::POWER2(FE::COS_K(this-> f,N),N);
    std::vector<cstm_float_t> cos2_coeff = c2->get_Taylor_coefficients(N);

    std::vector<cstm_float_t> T;
    T.push_back(c_tan(aux));

     T.push_back((f->get_coefficient(1))*c_pow(c_cos(aux),val(-2.0)));

             for (int k = 2; k < N ; ++k) {
            cstm_float_t termk = (f_coeffs[k]);//f->TAYLOR_K(k);
            cstm_float_t termj = val(0.0);
            for (int j = 1; j < k; ++j) {
                //aux =f->TAYLOR_K(k - j);
                termj += val(j) * T[j] * cos2_coeff[k-j];//Cos2.TAYLOR_K(k-j);

            }
            termj = val(-1.0)*termj / val(k);


            T.push_back(c_pow( c_cos(aux),val(-2.0) )*(termk+termj));
        }

        //------------------------------------------------
        this->coeffs = T;
}

cstm_float_t TAN_K::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);

        //if (k == 0){
        //    return c_tan(f->get_coefficient(0));
        //}
        //return coeffs[k];
        assert(k >= 0 && k < static_cast<int>(coeffs.size()));
        return coeffs[k];
}


ARCTAN_K::ARCTAN_K(std::shared_ptr<Function> f_, int N) : Function(N), f(f_) {
    //cstm_float_t aux = f->get_coefficient(0);
    std::vector<cstm_float_t> arg_f_coeff = f_->get_Taylor_coefficients(N);
    cstm_float_t aux = arg_f_coeff[0];


    auto P_function = FE::POWER2(f_,N);
    std::vector<cstm_float_t> pow_coeff = P_function ->get_Taylor_coefficients(N);

    std::vector<cstm_float_t> G, ARCT;
    G.push_back(val(1.0)+c_pow(aux,val(2.0)));

    //for(int i=1; i<=N; i++){
        //G.push_back( Pow.TAYLOR_K(i) );
    //    G.push_back( pow_coeff[i] );
    //}

    ARCT.push_back(c_atan(aux));
    cstm_float_t cte1 = val(1.0)/(val(1.0)+ c_pow(aux,val(2.0)));

        for (int k = 1; k < N ; ++k) {

            cstm_float_t termj = val(0.0);
            for (int j = 1; j < k; ++j) {
                //aux =f->TAYLOR_K(k - j);
                termj += val(j) * ARCT[j] * pow_coeff[k-j];

            }

            //ARCT.push_back( cte1 * (f->TAYLOR_K(k) - termj/val(k)));
            ARCT.push_back( cte1 * (arg_f_coeff[k] - termj/val(k)));

        }

        this->coeffs = ARCT;
}

cstm_float_t ARCTAN_K::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);
        if (k == 0){
            return c_atan(f->get_coefficient(0));
        }
        return coeffs[k];
}




POWER_K::POWER_K(cstm_float_t exponent_, std::shared_ptr<Function> f_, int N) : exponent(exponent_), Function(N), f(f_) {
    std::vector<cstm_float_t> f_Taylor = f_->get_Taylor_coefficients(N);
    cstm_float_t aux = f_Taylor[0];

    std::vector<cstm_float_t> P ;
    P.push_back(c_pow(aux,exponent_));

    for (int k = 1; k < N ; ++k) {

            cstm_float_t termj = val(0.0);
            for (int j = 0; j < k; ++j) {
                //aux =f->TAYLOR_K(k - j);
                termj += (exponent * val(k-j)-val(j) ) * P[j] * f_Taylor[k-j];

            }

            P.push_back( termj/(val(k)*aux) );

        }

        this->coeffs = P;

}

cstm_float_t POWER_K::get_coefficient(int k) const{
        //assert(k >= 0 && k < max_order);
        /*
        if (k == 0){
            return c_pow(f->get_coefficient(0),this->exponent);
        }
        return coeffs[k];
        */
        assert(k >= 0 && k < static_cast<int>(coeffs.size()));
    return coeffs[k];
}


SUM_K::SUM_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N) : Function(N), f(f_), g(g_) {
    std::vector<cstm_float_t> S;//, f_Taylor, g_Taylor;
    std::vector<cstm_float_t> f_Taylor = f->get_Taylor_coefficients(N);
    std::vector<cstm_float_t> g_Taylor = g->get_Taylor_coefficients(N);

    for (int k = 0; k < N ; ++k) {

    S.push_back( f_Taylor[k]+g_Taylor[k] );

    }

        this->coeffs = S;
}

cstm_float_t SUM_K::get_coefficient(int k) const{
    assert(k >= 0 && k < coeffs.size());
    return coeffs[k];
}



PROD_K::PROD_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N) : Function(N), f(f_), g(g_) {
    std::vector<cstm_float_t> P;//, f_Taylor, g_Taylor;
    std::vector<cstm_float_t> f_Taylor = f->get_Taylor_coefficients(N);
    std::vector<cstm_float_t> g_Taylor = g->get_Taylor_coefficients(N);


    this->N = N;

    for (int k = 0; k < N ; ++k) {
        cstm_float_t term_aux = val(0.0);
        for (int j = 0; j <= k ; ++j) {
            term_aux += f_Taylor[j]*g_Taylor[k-j];
        }
        P.push_back( term_aux );

    }

    this->coeffs = P;
}

cstm_float_t PROD_K::get_coefficient(int k) const{
    assert(k >= 0 && k < coeffs.size());
    return coeffs[k];
}


CTE_MULT::CTE_MULT(cstm_float_t ctev, std::shared_ptr<Function> f_, int N) : Function(N), cte_value(ctev), f(f_) {


    std::vector<cstm_float_t> S ;


    std::vector<cstm_float_t> arg_f_coeff = f -> get_Taylor_coefficients(N);

    for(int i=0; i<N ; i++){
       // S.push_back( f -> TAYLOR_K(i) *cte_value ); // f -> TAYLOR_K(i)
        S.push_back( arg_f_coeff[i] *ctev ); // f -> TAYLOR_K(i)
    }


        this->coeffs = S;
        this->N = N;
        this->cte_value = ctev;

}

cstm_float_t CTE_MULT::get_coefficient(int k) const{
    assert(k >= 0 && k < coeffs.size());
    return coeffs[k];
}

// -------------------------------------------------------------------------

cstm_float_t EvaluateTaylor(const std::shared_ptr<Function> &f_, cstm_float_t delta){
    cstm_float_t evaluation = val(0.0);
    int nn = f_->get_max_order();
    std::vector<cstm_float_t> tmp = f_->get_Taylor_coefficients(nn);
    //int nn = tmp.size();
    for (int j = 0; j<nn; j++ ) {
        evaluation += tmp[j]*c_pow(delta,val(j));
    }
    return evaluation;
}


std::vector<cstm_float_t> Get_FirstTerms(std::vector<cstm_float_t> series, cstm_float_t increment){
    std::vector<cstm_float_t> terms;
    cstm_float_t aux_beta0 = val(0.0);
    cstm_float_t aux_beta1 = val(0.0);
    cstm_float_t aux_beta2 = val(0.0);
    int N=series.size();
    for(int k=0; k<N;++k){
        aux_beta0 += series[k]*c_pow((increment),val(k*1.0));
        if(k!=0){
            aux_beta1 += val(k*1.0)*series[k]*c_pow((increment),val((k-1)*1.0));
        }
        if(k>1){
                aux_beta2 += val(k)*val(k-1)*series[k]*c_pow((increment),val(k-2));
            }
        //}
    }

    terms.push_back(aux_beta0);
    terms.push_back(aux_beta1);
    terms.push_back(aux_beta2);
    return terms;
}
// 
// MATH FUNCTIONS

template <typename T>
cstm_float_t c_pow(T value, T exponent) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return powq(value, exponent);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return pow(value,exponent);
        } else {
             return pow(value,exponent);
        }
    }
}

template cstm_float_t c_pow(float,float);
template cstm_float_t c_pow(double,double);
template cstm_float_t c_pow(__float128,__float128);



template <typename T>
cstm_float_t c_sin(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return sinq(value);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return sin(value);
        } else {
             return sin(value);
        }
    }
}

template cstm_float_t c_sin(float);
template cstm_float_t c_sin(double);
template cstm_float_t c_sin(__float128);




template <typename T>
cstm_float_t c_cos(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return cosq(value);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return cosq(value);
        } else {
             return cosq(value);
        }
    }
}

template cstm_float_t c_cos(float);
template cstm_float_t c_cos(double);
template cstm_float_t c_cos(__float128);


//
template <typename T>
cstm_float_t c_sqrt(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return sqrtq(value);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return sqrt(value);
        } else {
             return sqrt(value);
        }
    }
}

template cstm_float_t c_sqrt(float);
template cstm_float_t c_sqrt(double);
template cstm_float_t c_sqrt(__float128);


//
template <typename T>
cstm_float_t c_atan(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return atanq(value);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return atan(value);
        } else {
             return atan(value);
        }
    }
}

template cstm_float_t c_atan(float);
template cstm_float_t c_atan(double);
template cstm_float_t c_atan(__float128);



//
template <typename T>
cstm_float_t c_tan(T value) {
    if constexpr (std::is_same<T, __float128>::value ) {
        return tanq(value);

    } else {
        if constexpr (std::is_same<T, double>::value) {
             return tan(value);
        } else {
             return tan(value);
        }
    }
}

template cstm_float_t c_tan(float);
template cstm_float_t c_tan(double);
template cstm_float_t c_tan(__float128);
