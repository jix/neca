#ifndef NECA_GMP_HELPERS_H
#define NECA_GMP_HELPERS_H

#include <gmpxx.h>

static inline void set(mpz_class &rop, mpz_class const &op) {
    mpz_set(rop.get_mpz_t(), op.get_mpz_t());
}

static inline void neg(mpz_class &rop, mpz_class const &op) {
    mpz_neg(rop.get_mpz_t(), op.get_mpz_t());
}

static inline void neg(mpz_class &rop) {
    mpz_neg(rop.get_mpz_t(), rop.get_mpz_t());
}

static inline void sqrt(mpz_class &rop, mpz_class const &op) {
    mpz_sqrt(rop.get_mpz_t(), op.get_mpz_t());
}

static inline void sqrt(mpz_class &rop) {
    mpz_sqrt(rop.get_mpz_t(), rop.get_mpz_t());
}

static inline void set_ui(mpz_class &rop, unsigned long op) {
    mpz_set_ui(rop.get_mpz_t(), op);
}

static inline void set_si(mpz_class &rop, signed long op) {
    mpz_set_si(rop.get_mpz_t(), op);
}

static inline int set_str(mpz_class &rop, char const *str, int base) {
    return mpz_set_str(rop.get_mpz_t(), str, base);
}

static inline void mul(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_mul(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void mul(mpz_class &rop, mpz_class const &op2) {
    mpz_mul(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void mul_si(mpz_class &rop, mpz_class const &op1, long op2) {
    mpz_mul_si(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void mul_si(mpz_class &rop, long op2) {
    mpz_mul_si(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}


static inline void mul_2exp(mpz_class &rop, mpz_class const &op1, mp_bitcnt_t op2) {
    mpz_mul_2exp(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void mul_2exp(mpz_class &rop, mp_bitcnt_t op2) {
    mpz_mul_2exp(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline void add(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_add(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void add(mpz_class &rop, mpz_class const &op2) {
    mpz_add(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void add_ui(mpz_class &rop, mpz_class const &op1, unsigned long op2) {
    mpz_add_ui(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void add_ui(mpz_class &rop, unsigned long op2) {
    mpz_add_ui(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline void sub(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_sub(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void sub(mpz_class &rop, mpz_class const &op2) {
    mpz_sub(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void sub_ui(mpz_class &rop, mpz_class const &op1, unsigned long op2) {
    mpz_sub_ui(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void sub_ui(mpz_class &rop, unsigned long op2) {
    mpz_sub_ui(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline void mod(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_mod(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void mod(mpz_class &rop, mpz_class const &op2) {
    mpz_mod(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void invert(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_invert(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void invert(mpz_class &rop, mpz_class const &op2) {
    mpz_invert(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void addmul(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_addmul(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void submul(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_submul(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void divexact(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_divexact(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void divexact(mpz_class &rop, mpz_class const &op2) {
    mpz_divexact(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void tdiv_q_2exp(mpz_class &rop, mpz_class const &op1, mp_bitcnt_t op2) {
    mpz_tdiv_q_2exp(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void tdiv_q_2exp(mpz_class &rop, mp_bitcnt_t op2) {
    mpz_tdiv_q_2exp(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline void fdiv_q(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_fdiv_q(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void fdiv_q(mpz_class &rop, mpz_class const &op2) {
    mpz_fdiv_q(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void tdiv_q(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_tdiv_q(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void tdiv_q(mpz_class &rop, mpz_class const &op2) {
    mpz_tdiv_q(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline int cmp_ui(mpz_class const &op1, unsigned long op2) {
    return mpz_cmp_ui(op1.get_mpz_t(), op2);
}

static inline int cmp(mpz_class const &op1, mpz_class const &op2) {
    return mpz_cmp(op1.get_mpz_t(), op2.get_mpz_t());
}

static inline int cmpabs_ui(mpz_class const &op1, unsigned long op2) {
    return mpz_cmpabs_ui(op1.get_mpz_t(), op2);
}

static inline int cmpabs(mpz_class const &op1, mpz_class const &op2) {
    return mpz_cmpabs(op1.get_mpz_t(), op2.get_mpz_t());
}

static inline int sgn(mpz_class const &op) {
    return mpz_sgn(op.get_mpz_t());
}

static inline void powm(mpz_class &rop, mpz_class const &base, mpz_class const &exp, mpz_class const &mod) {
    mpz_powm(rop.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
}

static inline void powm_ui(mpz_class &rop, mpz_class const &base, unsigned long exp, mpz_class const &mod) {
    mpz_powm_ui(rop.get_mpz_t(), base.get_mpz_t(), exp, mod.get_mpz_t());
}

static inline void lcm(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_lcm(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void lcm(mpz_class &rop, mpz_class const &op2) {
    mpz_lcm(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void lcm_ui(mpz_class &rop, mpz_class const &op1, unsigned long op2) {
    mpz_lcm_ui(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline void lcm_ui(mpz_class &rop, unsigned long op2) {
    mpz_lcm_ui(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline void gcd(mpz_class &rop, mpz_class const &op1, mpz_class const &op2) {
    mpz_gcd(rop.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t());
}

static inline void gcd(mpz_class &rop, mpz_class const &op2) {
    mpz_gcd(rop.get_mpz_t(), rop.get_mpz_t(), op2.get_mpz_t());
}

static inline void gcdext(mpz_class &g, mpz_class &s, mpz_class &t, mpz_class const &a, mpz_class const &b) {
    mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}

static inline void gcdext(mpz_class &g, mpz_class &s, std::nullptr_t t, mpz_class const &a, mpz_class const &b) {
    mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t, a.get_mpz_t(), b.get_mpz_t());
}

static inline unsigned long gcd_ui(mpz_class &rop, mpz_class const &op1, unsigned long op2) {
    return mpz_gcd_ui(rop.get_mpz_t(), op1.get_mpz_t(), op2);
}

static inline unsigned long gcd_ui(mpz_class &rop, unsigned long op2) {
    return mpz_gcd_ui(rop.get_mpz_t(), rop.get_mpz_t(), op2);
}

static inline unsigned long ui_gcd_ui(mpz_class &rop, unsigned long op2) {
    return mpz_gcd_ui(nullptr, rop.get_mpz_t(), op2);
}

static inline mp_limb_t getlimbn(mpz_class const &op, mp_size_t n) {
    return mpz_getlimbn(op.get_mpz_t(), n);
}


#endif

