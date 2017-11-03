#ifndef NECA_ARITH_H
#define NECA_ARITH_H

#include <cstdint>
#include <gmpxx.h>

typedef std::uint64_t ulimb_t;
typedef std::int64_t limb_t;
typedef __uint128_t ulimb2_t;
typedef __int128_t limb2_t;
static constexpr unsigned limb_bits = 64;

template<unsigned W>
class bbn {
public:
    limb_t limbs[W];

    bool negative() const {
        bool result;
        for (unsigned i = 0; i < W; ++i) {
            result = (limbs[i] != 0) ? limbs[i] < 0 : result;
        }
        return result;
    }

    void set(mpz_class const &value) {
        unsigned length = mpz_size(value.get_mpz_t());
        ulimb_t const *value_limbs = mpz_limbs_read(value.get_mpz_t());
        if (length > W) {
            length = W;
        }
        ulimb_t invert = 0;
        limb_t carry = 0;

        if (mpz_sgn(value.get_mpz_t()) < 0) {
            invert = limb_t(-1);
            carry = 1;
        }

        for (unsigned i = 0; i < W; ++i) {
            ulimb_t value_limb = i < length ? value_limbs[i] : 0;
            value_limb ^= invert;

            limb2_t sum = (limb2_t)value_limb + (limb2_t)carry;

            limb_t low = sum;

            limbs[i] = low;

            sum -= low;

            carry = sum >> limb_bits;


        }
    }

    void into(mpz_class &value) const {
        ulimb_t *value_limbs = mpz_limbs_write(value.get_mpz_t(), W);

        ulimb_t invert = 0;
        limb_t carry = 0;

        bool is_negative = negative();

        if (is_negative) {
            invert = limb_t(-1);
            carry = -1;
        }

        int new_width = 0;

        for (unsigned i = 0; i < W; ++i) {
            ulimb_t value_limb;

            limb2_t sum = (limb2_t)limbs[i] + (limb2_t)carry;

            value_limb = limb_t(sum) ^ invert;

            value_limbs[i] = value_limb;

            new_width = value_limb != 0 ? i + 1 : new_width;

            carry = sum >> limb_bits;
        }

        if (is_negative) {
            new_width = -new_width;
        }

        mpz_limbs_finish(value.get_mpz_t(), new_width);
    }

    double get_d() const {
        double r = 0;
        for (unsigned i = 0; i < W; ++i) {
            r *= 5.421010862427522e-20;
            r += limbs[i];
        }
        return r;
    }

    static void get_d_scaled4(
            double &r0, double &r1, double &r2, double &r3,
            bbn const &v0, bbn const &v1, bbn const &v2, bbn const &v3) {

        r0 = 0;
        r1 = 0;
        r2 = 0;
        r3 = 0;

        double scale = 5.421010862427522e-20;

        for (unsigned i = 0; i < W; ++i) {

            limb_t l0 = v0.limbs[i];
            limb_t l1 = v1.limbs[i];
            limb_t l2 = v2.limbs[i];
            limb_t l3 = v3.limbs[i];

            if (l0 || l1 || l2 || l3) {
                r0 *= scale;
                r1 *= scale;
                r2 *= scale;
                r3 *= scale;


                r0 += l0;
                r1 += l1;
                r2 += l2;
                r3 += l3;


                scale = 5.421010862427522e-20;
            } else {
                scale *= 5.421010862427522e-20;
            }
        }
    }

    __attribute__((noinline))
    static void apply_mat4(
            bbn &v00, bbn &v01, bbn &v10, bbn &v11,
            limb_t m00, limb_t m01, limb_t m10, limb_t m11) {

        // v00 = m00 * v00 + m01 * v10
        // v01 = m00 * v01 + m01 * v11
        // v10 = m10 * v00 + m11 * v10
        // v11 = m10 * v01 + m11 * v11

        limb_t cx_00 = 0, cx_01 = 0, cx_10 = 0, cx_11 = 0;
        limb_t cy_00 = 0, cy_01 = 0, cy_10 = 0, cy_11 = 0;
        limb_t cz_00 = 0, cz_01 = 0, cz_10 = 0, cz_11 = 0;

        for (unsigned i = 0; i < W; ++i) {
            limb2_t a00 = (limb2_t)m00 * (limb2_t)v00.limbs[i];
            limb2_t b00 = (limb2_t)m01 * (limb2_t)v10.limbs[i];

            limb2_t a01 = (limb2_t)m00 * (limb2_t)v01.limbs[i];
            limb2_t b01 = (limb2_t)m01 * (limb2_t)v11.limbs[i];


            limb2_t a10 = (limb2_t)m10 * (limb2_t)v00.limbs[i];
            limb2_t b10 = (limb2_t)m11 * (limb2_t)v10.limbs[i];

            limb2_t a11 = (limb2_t)m10 * (limb2_t)v01.limbs[i];
            limb2_t b11 = (limb2_t)m11 * (limb2_t)v11.limbs[i];

            limb_t al_00 = a00; a00 -= al_00;
            limb_t al_01 = a01; a01 -= al_01;
            limb_t al_10 = a10; a10 -= al_10;
            limb_t al_11 = a11; a11 -= al_11;

            limb_t bl_00 = b00; b00 -= bl_00;
            limb_t bl_01 = b01; b01 -= bl_01;
            limb_t bl_10 = b10; b10 -= bl_10;
            limb_t bl_11 = b11; b11 -= bl_11;

            limb2_t r00 = (limb2_t)al_00 + (limb2_t)bl_00 + (limb2_t)cx_00 + (limb2_t)cy_00 + (limb2_t)cz_00;
            limb2_t r01 = (limb2_t)al_01 + (limb2_t)bl_01 + (limb2_t)cx_01 + (limb2_t)cy_01 + (limb2_t)cz_01;
            limb2_t r10 = (limb2_t)al_10 + (limb2_t)bl_10 + (limb2_t)cx_10 + (limb2_t)cy_10 + (limb2_t)cz_10;
            limb2_t r11 = (limb2_t)al_11 + (limb2_t)bl_11 + (limb2_t)cx_11 + (limb2_t)cy_11 + (limb2_t)cz_11;


            limb_t rl_00 = r00; r00 -= rl_00;
            limb_t rl_01 = r01; r01 -= rl_01;
            limb_t rl_10 = r10; r10 -= rl_10;
            limb_t rl_11 = r11; r11 -= rl_11;

            v00.limbs[i] = rl_00;
            v01.limbs[i] = rl_01;
            v10.limbs[i] = rl_10;
            v11.limbs[i] = rl_11;

            cx_00 = r00 >> limb_bits; cy_00 = a00 >> limb_bits; cz_00 = b00 >> limb_bits;
            cx_01 = r01 >> limb_bits; cy_01 = a01 >> limb_bits; cz_01 = b01 >> limb_bits;
            cx_10 = r10 >> limb_bits; cy_10 = a10 >> limb_bits; cz_10 = b10 >> limb_bits;
            cx_11 = r11 >> limb_bits; cy_11 = a11 >> limb_bits; cz_11 = b11 >> limb_bits;
        }
    }

};

/*

inline void limb_mul(
        limb_t &high, limb_t &low, limb_t a, limb_t b) {
    __uint128_t wide = (__uint128_t)a * (__uint128_t)b;
    high = wide >> limb_bits;
    low = wide;
}

inline void limb_mul_add(
        limb_t &high, limb_t &low, limb_t a, limb_t b, limb_t c) {
    __uint128_t wide = (__uint128_t)a * (__uint128_t)b;
    limb_t h = wide >> limb_bits;
    limb_t l = wide;

    l += c;
    h += l < c;

    low = l;
    high = h;
}

inline void limb_smul(
        limb_t &high, limb_t &low, limb_t a, limb_t b) {
    __uint128_t wide = (__uint128_t)(
        (__int128_t)(slimb_t)a * (__int128_t)(slimb_t)b);
    high = wide >> limb_bits;
    low = wide;
}

inline void limb_smul_add(
        limb_t &high, limb_t &low, limb_t a, limb_t b, limb_t c) {
    __uint128_t wide = (__uint128_t)(
        (__int128_t)(slimb_t)a * (__int128_t)(slimb_t)b);
    limb_t h = wide >> limb_bits;
    limb_t l = wide;

    l += c;
    h += l < c;

    low = l;
    high = h;
}


inline void limb_add(
        limb_t &high, limb_t &low, limb_t a, limb_t b) {
    low = a + b;
    high = low < a;
}

inline void limb_addc(
        limb_t &high, limb_t &low, limb_t a, limb_t b, limb_t c) {
    limb_t l = a + b;
    limb_t h = l < a;
    l += c;
    h += l < c;
    low = l;
    high = h;
}

inline void limb_sub(
        limb_t &high, limb_t &low, limb_t a, limb_t b) {
    low = a - b;
    high = low > a;
}

inline void limb_subb(
        limb_t &high, limb_t &low, limb_t a, limb_t b, limb_t c) {
    limb_t l = a - b;
    limb_t h = l > a;
    limb_t lpb = l;
    l -= c;
    h += l > lpb;
    low = l;
    high = h;
}

inline limb_t limb_signext(limb_t a)
{
    return ((slimb_t)a) >> (limb_bits - 1);
}

*/


#endif
