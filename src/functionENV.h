#ifndef FUNCTIONENV_H
#define FUNCTIONENV_H

#include <memory>
#include "PowerSeries.h"

namespace FE {

    inline std::shared_ptr<Function> Variable(cstm_float_t p, int N) {
        return std::make_shared<::Variable>(p, N);
    }

    inline std::shared_ptr<Function> CTE(cstm_float_t value, int N) {
        return std::make_shared<::CTE>(value, N);
    }

    inline std::shared_ptr<Function> POWER2(std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::POWER2>(f_, N);
    }

    inline std::shared_ptr<Function> SIN_K(std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::SIN_K>(f_, N);
    }

    inline std::shared_ptr<Function> COS_K(std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::COS_K>(f_, N);
    }

    inline std::shared_ptr<Function> TAN_K(std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::TAN_K>(f_, N);
    }

    inline std::shared_ptr<Function> ARCTAN_K(std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::ARCTAN_K>(f_, N);
    }

    inline std::shared_ptr<Function> POWER_K(cstm_float_t exponent_, std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::POWER_K>(exponent_, f_, N);
    }


    inline std::shared_ptr<Function> SUM_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N) {
        return std::make_shared<::SUM_K>(f_, g_, N);
    }

    inline std::shared_ptr<Function> PROD_K(std::shared_ptr<Function> f_, std::shared_ptr<Function> g_, int N) {
        return std::make_shared<::PROD_K>(f_, g_, N);
    }

    inline std::shared_ptr<Function> CTE_MULT(cstm_float_t ctev, std::shared_ptr<Function> f_, int N) {
        return std::make_shared<::CTE_MULT>(ctev, f_, N);
    }
}

#endif
