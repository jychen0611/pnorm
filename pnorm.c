#include<stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define DIMENSION 3
#define TEST_CASE1 {1.0, 2.0, 3.0}
#define TEST_CASE2 {1.0, -2.0, -3.0}
#define TEST_CASE3 {0.0, -1.0, -7.0}

#define iswap(x, y) ((x) ^= (y), (y) ^= (x), (x) ^= (y))
static inline int64_t getbit(int64_t value, int n) { return (value >> n) & 1; }
static int32_t imul24(int32_t a, int32_t b);
static int32_t idiv24(int32_t a, int32_t b);

/* count_leading_zeros */
uint32_t count_leading_zeros(uint32_t x) {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    /* count ones (population count) */
    x -= ((x >> 1) & 0x55555555);
    x = ((x >> 2) & 0x33333333) + (x & 0x33333333);
    x = ((x >> 4) + x) & 0x0f0f0f0f;
    x += (x >> 8);
    x += (x >> 16);

    return (32 - (x & 0x7f));
}
/* my_clz */
static inline int my_clz(uint32_t x) {
    int count = 0;
    for (int i = 31; i >= 0; --i) {
        if (x & (1U << i))
            break;
        count++;
    }
    return count;
}

/* fmul32 */
float fmul32(float a, float b)
{
    int32_t ia = *(int32_t *)&a, ib = *(int32_t *)&b;
    if (ia == 0 || ib == 0) return 0;
    
    //int32_t s_a = (ia & 0x80000000) >> 31;
    //int32_t s_b = (ib & 0x80000000) >> 31;
    int32_t exp_a = (ia & 0x7F800000) >> 23;
    int32_t exp_b = (ib & 0x7F800000) >> 23;
    int32_t fra_a = ia & 0x7FFFFF;
    int32_t fra_b = ib & 0x7FFFFF;
    /* TODO: Special values like NaN and INF */
    if((exp_a == 255 && fra_a == 0 ) || (exp_b == 255 && fra_b == 0) ){
        printf("Infinity value!\n");
    }
    if((exp_a == 255 && fra_a != 0 ) || (exp_b == 255 && fra_b != 0) ){
        printf("NaN value!\n");
    }

    /* mantissa */
    int32_t ma = (ia & 0x7FFFFF) | 0x800000;
    int32_t mb = (ib & 0x7FFFFF) | 0x800000;

    int32_t sea = ia & 0xFF800000;
    int32_t seb = ib & 0xFF800000;

    /* result of mantissa */
    int32_t m = imul24(ma, mb);
    int32_t mshift = getbit(m, 24);
    m >>= mshift;

    //int32_t r = (~(s_a ^ s_b) << 31) + ((exp_a + (exp_b - 0x7F)) << 23) + m;
    int32_t r = ((sea - 0x3f800000 + seb) & 0xFF800000) + m - (0x800000 & -!mshift);
    int32_t ovfl = (r ^ ia ^ ib) >> 31;
    r = r ^ ( (r ^ 0x7f800000) & ovfl);
    return *(float *)&r;
}

/* imul24 */
static int32_t imul24(int32_t a, int32_t b)
{
    uint32_t r = 0;
    for (; b; b >>= 4){
        r = (r >> 1) + (a & -getbit(b, 0));
        r = (r >> 1) + (a & -getbit(b, 1));
        r = (r >> 1) + (a & -getbit(b, 2));
        r = (r >> 1) + (a & -getbit(b, 3));
    }
    return r;
}

/* fdiv32 */
float fdiv32(float a, float b)
{
    int32_t ia = *(int32_t *)&a, ib = *(int32_t *)&b;
    if (a == 0) return a;
    if (b == 0) return *(float*)&(int){0x7f800000};
    /* mantissa */
    int32_t ma = (ia & 0x7FFFFF) | 0x800000;
    int32_t mb = (ib & 0x7FFFFF) | 0x800000;

    /* sign and exponent */
    int32_t sea = ia & 0xFF800000;
    int32_t seb = ib & 0xFF800000;

    /* result of mantissa */
    int32_t m = idiv24(ma, mb);
    int32_t mshift = !getbit(m, 31);
    m <<= mshift;

    int32_t r = ((sea - seb + 0x3f800000) - (0x800000 & -mshift)) | (m & 0x7fffff00) >> 8;
    int32_t ovfl = (ia ^ ib ^ r) >> 31;
    r = r ^ ((r ^ 0x7f800000) & ovfl);
    
    return *(float *) &r;
    // return a / b;
}

/* idiv24 */
static int32_t idiv24(int32_t a, int32_t b) {
    uint32_t r = 0;
    for (int i = 0; i < 32; i++) {
        a -= b;
        r = (r << 1) | (a >= 0);
        a = (a + (b & -(a < 0))) << 1;
    }

    return r;
}

/* fadd32 */
float fadd32(float a, float b) {
    int32_t ia = *(int32_t *)&a, ib = *(int32_t *)&b;

    int32_t cmp_a = ia & 0x7fffffff;
    int32_t cmp_b = ib & 0x7fffffff;

    if (cmp_a < cmp_b)
        iswap(ia, ib);
    /* exponent */
    int32_t ea = (ia >> 23) & 0xff;
    int32_t eb = (ib >> 23) & 0xff;

    /* mantissa */
    int32_t ma = (ia & 0x7fffff) | 0x800000;
    int32_t mb = (ib & 0x7fffff) | 0x800000;

    int32_t align = (ea - eb > 24) ? 24 : (ea - eb);

    mb >>= align;
    if ((ia ^ ib) >> 31) {
        ma -= mb;
    } else {
        ma += mb;
    }

    int32_t clz = my_clz(ma);
    int32_t shift = 0;
    if (clz <= 8) {
        shift = 8 - clz;
        ma >>= shift;
        ea += shift;
    } else {
        shift = clz - 8;
        ma <<= shift;
        ea -= shift;
    }

    int32_t r = (ia & 0x80000000) | ea << 23 | (ma & 0x7fffff);
    //float tr = a + b;
    return *(float *)&r;
}

/* f2i32 */
int f2i32(int x) {
    int32_t a = *(int *)&x;
    int32_t ma = (a & 0x7FFFFF) | 0x800000;
    int32_t ea = ((a >> 23) & 0xFF) - 127;
    if (ea < 0)
        return 0;
    else if (ea <= 23)
        ma >>= (23 - ea);
    else
        ma <<= (ea - 23);

    return ma;
}

/* i2f32 */
int i2f32(int x) {
    if (x == 0) return 0;

    int32_t s = x & 0x80000000;
    if (s) x = -x;

    int32_t clz = my_clz(x);
    int32_t e = 31 - clz + 127;

    if (clz <= 8) {
        x >>= 8 - clz;
    } else {
        x <<= clz - 8;
    }

    int r = s | e << 23 | (x & 0x7fffff);
    return r;
}

/* fabsf */
static inline float fabsf(float x){
    uint32_t i = *(uint32_t *)&x;
    i &= 0x7FFFFFFF;
    x = *(float *)&i;
    return x;
}

/* my_exp */
static inline float my_exp(float x) {
    float result = 1.0f; 
    float term = 1.0f;   
    for(int n = 1; n <= 5; n++){ 
        term = fmul32(term, fdiv32(x, (float)n));
        result = fadd32(result, term); 
    }
    return result;
}

/* my_ln */
static inline float my_ln(float x) {
    if(x <= 0.0f)
        return -1.0f;
    if(x == 1.0f)
        return 0.0f;
    if(x <= 2.8f)
        return 1.0f;

    float result = 0.0f;
    float term = fdiv32(fadd32(x, -1.0f), fadd32(x, 1.0f));
    float term_squared = fmul32(term, term);
    float current_term = term;
    for(int n = 1; n <= 100; n += 2){
        result = fadd32(result, fdiv32(current_term, n)); 
        current_term = fmul32(current_term, term_squared);
    }

    return fmul32(2.0f, result);
}

/* my_pow */
static inline float my_pow(float base, float exponent){
    if(base == 0.0f && exponent < 0.0f)
        return 0.0f;

    if(exponent == 0.0f)
        return 1.0f;

    float result = 1.0f;
    int exp_int = (int)exponent;
    float fraction = fadd32(exponent, (float)-exp_int);

    if(exp_int < 0){
        base = fdiv32(1.0f, base);
        exp_int = -exp_int;
    }

    for(int i = 0; i < exp_int; ++i)
        result = fmul32(result, base);

    if(fraction > 0)
        result = fmul32(result, my_exp(fmul32(fraction, my_ln(base))));

    return result;
}

/* my_root */
static inline float my_root(float x, int n){
    return my_pow(x, fdiv32(1.0f, (float)n));
}

/* find_norm */
static inline float find_norm(float* v, int n, int p){
    float res = 0.0f;
    for(int i=0;i<n;++i){
        float pow_val =  my_pow(fabsf(v[i]), p);
        res = fadd32(res, pow_val);
    }
    res = my_root(res, p);
    return res;
}

/* main */
int main(void){
    float vec1[] = TEST_CASE1;
    printf("TEST_CASE1 {1.0, 2.0, 3.0}\n");
    for(int p=1;p<=2;++p){
        float pnorm = find_norm(vec1, DIMENSION, p);
        printf("%d-norm : %f\n", p, pnorm);
    }
    float vec2[] = TEST_CASE2;
    printf("TEST_CASE2 {1.0, -2.0, 3.0}\n");
    for(int p=1;p<=2;++p){
        float pnorm = find_norm(vec2, DIMENSION, p);
        printf("%d-norm : %f\n", p, pnorm);
    }
    float vec3[] = TEST_CASE3;
    printf("TEST_CASE3 {0.0, -1.0, -7.0}\n");
    for(int p=1;p<=2;++p){
        float pnorm = find_norm(vec3, DIMENSION, p);
        printf("%d-norm : %f\n", p, pnorm);
    }
    
    return 0;
}

