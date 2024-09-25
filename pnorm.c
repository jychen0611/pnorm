#include<stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define VECTOR {1.0, 2.0, 3.0}
#define DIMENSION 3

static inline float fabsf(float x){
    uint32_t i = *(uint32_t *)&x;
    i &= 0x7FFFFFFF;
    x = *(float *)&i;
    return x;
}

static inline float my_exp(float x) {
    float result = 1.0f; 
    float term = 1.0f;   
    for(int n = 1; n <= 20; n++){ 
        term *= x / n; 
        result += term; 
    }
    return result;
}

static inline float my_ln(float x) {
    if(x <= 0.0f)
        return -1.0f;

    float result = 0.0f;
    float term = (x - 1) / (x + 1);
    float term_squared = term * term;
    float current_term = term;
    for(int n = 1; n <= 100; n += 2){
        result += current_term / n;
        current_term *= term_squared;
    }

    return 2 * result;
}

static inline float my_pow(float base, float exponent){
    if(base == 0.0f && exponent < 0.0f)
        return 0.0f;

    if(exponent == 0.0f)
        return 1.0f;

    float result = 1.0f;
    int exp_int = (int)exponent;
    float fraction = exponent - exp_int;

    if(exp_int < 0){
        base = 1.0f / base;
        exp_int = -exp_int;
    }

    for(int i = 0; i < exp_int; ++i)
        result *= base;

    if(fraction > 0)
        result *= my_exp(fraction * my_ln(base));

    return result;
}

static inline float my_root(float x, int n){
    return my_pow(x, 1.0 / n);
}

static inline float find_norm(float* v, int n, int p){
    float res = 0.0f;
    for(int i=0;i<n;++i){
        res +=  my_pow(fabsf(v[i]), p);
    }
    res = my_root(res, p);
    return res;
}

int main(void){
    float vec[] = VECTOR;
    for(int p=1;p<=2;++p){
        float pnorm = find_norm(vec, DIMENSION, p);
        printf("%d-norm : %f\n", p, pnorm);
    }

    return 0;
}

