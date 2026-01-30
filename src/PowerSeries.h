/*
Posc++ - Library for Power-Series Composition
The files in this project are maintained in the GitHub repository Posc++, available at
https://github.com/V3du4rd0/posc-.
*/

#ifndef PowerSeries_H
#define PowerSeries_H

#include "config.h"
#include <quadmath.h>

#include <type_traits>

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <map>
#include <sys/stat.h>
#include <limits>
#include <math.h>
#include <memory>
#include <cassert>


class MathExpression;
class Function;
class Constant;
class Variable;


#ifdef USE_FLOAT128
    typedef __float128 cstm_float_t;
#elif defined USE_DOUBLE
    typedef double cstm_float_t;
#else
    typedef float cstm_float_t;
#endif


extern const cstm_float_t Pi;


template <typename T>
void print_result(T value);


cstm_float_t val(double value);


//

class Function {
protected:
    int max_order;

public:
    Function(int N);
    virtual ~Function();

    virtual cstm_float_t get_coefficient(int k) const = 0;
    virtual std::vector<cstm_float_t> get_Taylor_coefficients(int N) const;
    int get_max_order() const { return max_order; }
};


class Variable : public Function {
    cstm_float_t p;
    std::vector<cstm_float_t> coeffs;

public:
    Variable(cstm_float_t p_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class CTE : public Function {
    cstm_float_t value;
    std::vector<cstm_float_t> coeffs;
public:
    CTE(cstm_float_t val, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class POWER2 : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
public:
    POWER2(std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};




class SIN_K : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
public:
    SIN_K(std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class COS_K : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
public:
    COS_K(std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class TAN_K : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
public:
    TAN_K(std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class ARCTAN_K : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
public:
    ARCTAN_K(std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class POWER_K : public Function {
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
    cstm_float_t exponent;
public:
    POWER_K(cstm_float_t exponent_, std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class SUM_K : public Function{
    std::shared_ptr<Function> f;
    std::shared_ptr<Function> g;
    std::vector<cstm_float_t> coeffs;

public:
    SUM_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class PROD_K : public Function{
    std::shared_ptr<Function> f;
    std::shared_ptr<Function> g;
    std::vector<cstm_float_t> coeffs;
    int N;

public:
    PROD_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


class CTE_MULT : public Function{
    std::shared_ptr<Function> f;
    std::vector<cstm_float_t> coeffs;
    cstm_float_t cte_value;
    int N;

public:
    CTE_MULT(cstm_float_t ctev, std::shared_ptr<Function> f_, int N);
    cstm_float_t get_coefficient(int k) const override;
};


//
cstm_float_t EvaluateTaylor(const std::shared_ptr<Function> &f_, cstm_float_t delta);

std::vector<cstm_float_t> Get_FirstTerms(std::vector<cstm_float_t> series, cstm_float_t increment);
// 

template <typename T>
cstm_float_t c_pow(T value, T exponent);

template <typename T>
cstm_float_t c_sin(T value);

template <typename T>
cstm_float_t c_cos(T value);

template <typename T>
cstm_float_t c_sqrt(T value);

template <typename T>
cstm_float_t c_atan(T value);

template <typename T>
cstm_float_t c_tan(T value);

#endif
