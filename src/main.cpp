#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "gmp_helpers.h"
#include "arith.h"

using std::cout;
using std::endl;

template<class T>
using vec = std::array<T, 2>;

template<class T>
using mat = vec<vec<T>>;

#if defined(__GNUC__)
#define perf_noinline __attribute__((noinline))
#else
#define perf_noinline
#endif

class neca {
public:
    mpz_class n;
    mpz_class gen;

    mpz_class s;
    mpz_class half_s;
    mpz_class n_inv;
    mpz_class gen_inv;
    mpz_class dlog_n;
    mpz_class order;

    mpz_class sqrt_n_div;
    mpz_class max_step_sqnorm;

    mpz_class p, q;
    mpz_class p_mod, q_mod;

    mpz_class eval;

    mpz_class line_min, line_range;
    mpz_class line_sign_min, line_sign_range;

    mat<mpz_class> basis;
    vec<mpz_class> origin;

    mpz_class tmp[4];

    neca(mpz_class const &n, mpz_class const &gen,
            std::vector<unsigned long> const &primes) {
        this->n = n;
        this->gen = gen;

        s = 1;
        mpz_class x, y, z;
        dlog_n = 0;
        order = 1;

        for (unsigned long p_i : primes) {
            s *= p_i;

            x = gen % p_i;
            y = n % p_i;
            unsigned order_p;
            unsigned dlog_p = 0;
            for (order_p = 1; x != 1; ++order_p) {
                if (x == y) {
                    dlog_p = order_p;
                }
                x *= gen;
                x %= p_i;
            }
            z = order_p;
            gcdext(x, z, y, z, order);

            dlog_n *= z * order_p;
            dlog_n += dlog_p * y * order;

            order *= order_p;
            dlog_n /= x;
            order /= x;
            mod(dlog_n, order);
        }

        powm(x, gen, dlog_n, s);

        x -= n;
        mod(x, s);

        if (x != 0) {
            order = 0;
            dlog_n = -1;
        }

        half_s = s >> 1;

        sqrt_n_div = sqrt(n);
        sqrt_n_div /= s;

        max_step_sqnorm = sqrt_n_div;
        max_step_sqnorm *= max_step_sqnorm;

        invert(n_inv, n, s);
        invert(gen_inv, gen, s);
    }

    perf_noinline
    void jump_guess(mpz_class const &i) {
        powm(p_mod, gen, i, s);
        powm(q_mod, gen_inv, i, s);
        mul(q_mod, n);
        mod(q_mod, s);
    }


    perf_noinline
    void next_guess() {
        mul(p_mod, gen);
        mul(q_mod, gen_inv);
        mod(p_mod, s);
        mod(q_mod, s);
    }

    perf_noinline
    void build_lattice() {
        mpz_class &q_inv = tmp[0];

        mul(q_inv, p_mod, n_inv);

        set(origin[0], n);
        submul(origin[0], p_mod, q_mod);
        divexact(origin[0], s);
        mul(origin[0], q_inv);
        mod(origin[0], s);

        set_ui(origin[1], 0);

        mul(basis[0][0], q_inv, p_mod);
        neg(basis[0][0]);
        mod(basis[0][0], s);

        set_ui(basis[0][1], 1);
        set(basis[1][0], s);
        set_ui(basis[1][1], 0);

        reduce_lattice();

        sub(tmp[0], sqrt_n_div, origin[0]);
        sub(tmp[1], sqrt_n_div, origin[1]);

        approx_cvp_inverse(tmp[2], tmp[3], tmp[0], tmp[1]);
        lattice_vec(tmp[0], tmp[1], tmp[2], tmp[3]);

        add(origin[0], tmp[0]);
        add(origin[1], tmp[1]);
    }

    perf_noinline
    bool search_factors() {
        build_lattice();

        mul(tmp[0], basis[0][0], basis[0][0]);
        addmul(tmp[0], basis[0][1], basis[0][1]);

        mul(tmp[1], basis[1][0], basis[1][0]);
        addmul(tmp[1], basis[1][1], basis[1][1]);

        vec<bool> search;

        search[0] = cmpabs(tmp[0], max_step_sqnorm) <= 0;
        search[1] = cmpabs(tmp[1], max_step_sqnorm) <= 0;

        if (!search[0] && !search[1]) {
            set(p, p_mod);
            addmul(p, s, origin[0]);
            set(q, q_mod);
            addmul(q, s, origin[1]);

            return eval_pq() == 0;
        } else if (search[0] && search[1]) {
            #pragma omp critical(user_io)
            {
                cout << "\nThis should not happen for 512-bit keys!\n" <<
                    n << "\n" <<
                    p_mod << endl;

                cout << (max_step_sqnorm / tmp[0]) << endl;
                cout << (max_step_sqnorm / tmp[1]) << endl;

                exit(1);
            }
        } else {
            if (!search[0]) {
                std::swap(basis[0][0], basis[1][0]);
                std::swap(basis[0][1], basis[1][1]);
                std::swap(tmp[0], tmp[1]);
            }

            tdiv_q(line_min, max_step_sqnorm, tmp[0]);
            sqrt(line_min);
            mul_2exp(line_range, line_min, 1);
            neg(line_min);

            mul(basis[0][0], s);
            mul(basis[0][1], s);

            mul(origin[0], s);
            mul(origin[1], s);

            add(origin[0], p_mod);
            add(origin[1], q_mod);

            return search_origin_line();
        }
    }

    perf_noinline
    bool search_origin_line() {
        mpz_class &mid = tmp[0];
        mpz_class &half = tmp[1];

        line_set_pq(line_min);

        int sign_min = eval_pq();
        if (sign_min == 0) {
            return true;
        }

        int dir_min = eval_d_pq();
        if (dir_min == 0) {
            return false;
        }

        line_move_pq(line_range);

        int sign_max = eval_pq();
        if (sign_max == 0) {
            return true;
        }

        int dir_max = eval_d_pq();
        if (dir_max == 0) {
            return false;
        }

        if (sign_min != sign_max) {
            if (sign_min > sign_max) {
                add(line_sign_min, line_min, line_range);
                neg(line_sign_range, line_range);
            } else {
                set(line_sign_min, line_min);
                set(line_sign_range, line_range);
            }
            return search_origin_line_sign();
        } else if (dir_min == dir_max) {
            return false;
        }

        int sign_outside = sign_min;

        if (dir_min > dir_max) {
            add(line_min, line_range);
            neg(line_range);
        }

        while (true) {
            tdiv_q_2exp(half, line_range, 1);
            if (sgn(half) == 0) {
                return false;
            }

            add(mid, line_min, half);

            line_set_pq(mid);

            int sign_mid = eval_pq();
            if (sign_mid == 0) {
                return true;
            }

            if (sign_mid != sign_outside) {
                if (sign_mid > sign_outside) {
                    set(line_sign_min, line_min);
                    set(line_sign_range, half);
                    if (search_origin_line_sign()) {
                        return true;
                    }
                    add(line_sign_min, line_min, line_range);
                    sub(line_sign_range, half, line_range);
                    return search_origin_line_sign();
                } else {
                    set(line_sign_min, mid);
                    neg(line_sign_range, half);
                    if (search_origin_line_sign()) {
                        return true;
                    }
                    set(line_sign_min, mid);
                    sub(line_sign_range, line_range, half);
                    return search_origin_line_sign();
                }
            }

            int dir_mid = eval_d_pq();
            if (dir_mid == 0) {
                return false;
            }

            if (dir_mid > 0) {
                set(line_range, half);
            } else {
                sub(line_range, half);
                set(line_min, mid);
            }
        }
    }

    perf_noinline
    bool search_origin_line_sign() {
        mpz_class &mid = tmp[0];
        mpz_class &half = tmp[1];

        while (true) {
            tdiv_q_2exp(half, line_sign_range, 1);
            if (sgn(half) == 0) {
                return false;
            }

            add(mid, line_sign_min, half);

            line_set_pq(mid);

            int sign_mid = eval_pq();
            if (sign_mid == 0) {
                return true;
            }

            if (sign_mid > 0) {
                set(line_sign_range, half);
            } else {
                sub(line_sign_range, half);
                set(line_sign_min, mid);
            }
        }
    }

    void line_set_pq(mpz_class const &i) {
        mul(p, basis[0][0], i);
        mul(q, basis[0][1], i);
        add(p, origin[0]);
        add(q, origin[1]);
    }

    void line_move_pq(mpz_class const &i) {
        addmul(p, basis[0][0], i);
        addmul(q, basis[0][1], i);
    }

    int eval_pq() {
        mul(eval, p, q);
        int diff = cmp(eval, n);
        return (0 < diff) - (diff < 0);
    }

    int eval_d_pq(int axis = 0) {
        mul(eval, p, basis[axis][1]);
        addmul(eval, q, basis[axis][0]);
        return sgn(eval);
    }

    perf_noinline
    void reduce_lattice() {
        typedef bbn<3> coord_t;

        mat<coord_t> bb_basis;

        bb_basis[0][0].set(basis[0][0]);
        bb_basis[0][1].set(basis[0][1]);
        bb_basis[1][0].set(basis[1][0]);
        bb_basis[1][1].set(basis[1][1]);

        bool done = false;
        for (unsigned counter = 0; !done; ++counter) {
            if (counter > 30) {
                fallback_reduce_lattice();
                return;
            }

            mat<double> db;
            mat<limb_t> trans = {
                {{1, 0}, {0, 1}}
            };

            coord_t::get_d_scaled4(
                db[0][0],
                db[0][1],
                db[1][0],
                db[1][1],
                bb_basis[0][0],
                bb_basis[0][1],
                bb_basis[1][0],
                bb_basis[1][1]
            );

            vec<double> lengths;

            lengths[0] = db[0][0] * db[0][0] + db[0][1] * db[0][1];
            lengths[1] = db[1][0] * db[1][0] + db[1][1] * db[1][1];

            double trans_size = 1.0;

            for (unsigned subcounter = 0; subcounter < 100; ++subcounter) {
                double inner = db[0][0] * db[1][0] + db[0][1] * db[1][1];

                int i = lengths[1] < lengths[0];
                double fa = inner / lengths[i];
                double ia = std::floor(fa + 0.5);
                trans_size *= std::abs(ia);

                limb_t a = ia;
                if (a != ia) {
                    fallback_reduce_lattice();
                    return;
                }

                if (subcounter > 0 && trans_size > double(1l << 29)) {
                    break;
                }

                if (a == 0) {
                    done = true;
                    break;
                }

                db[i ^ 1][0] -= db[i ^ 0][0] * ia;
                db[i ^ 1][1] -= db[i ^ 0][1] * ia;

                lengths[i ^ 1] =
                    db[i ^ 1][0] * db[i ^ 1][0] + db[i ^ 1][1] * db[i ^ 1][1];


                trans[i ^ 1][0] -= trans[i ^ 0][0] * a;
                trans[i ^ 1][1] -= trans[i ^ 0][1] * a;
            }

            coord_t::apply_mat4(
                bb_basis[0][0], bb_basis[0][1],
                bb_basis[1][0], bb_basis[1][1],
                trans[0][0], trans[0][1],
                trans[1][0], trans[1][1]
            );
        }

        bb_basis[0][0].into(basis[0][0]);
        bb_basis[0][1].into(basis[0][1]);
        bb_basis[1][0].into(basis[1][0]);
        bb_basis[1][1].into(basis[1][1]);
    }

    perf_noinline
    void fallback_reduce_lattice() {
        mpz_class &a = tmp[0];
        mpz_class &b = tmp[1];

        while (true) {
            for (int j = 0; j < 2; ++j) {
                int i = j ^ 0;

                mul(a, basis[i ^ 0][0], basis[i ^ 1][0]);
                addmul(a, basis[i ^ 0][1], basis[i ^ 1][1]);

                mul(b, basis[i ^ 0][0], basis[i ^ 0][0]);
                addmul(b, basis[i ^ 0][1], basis[i ^ 0][1]);

                mul_2exp(a, 1);
                add(a, b);
                mul_2exp(b, 1);

                fdiv_q(a, b);

                if (cmpabs_ui(a, 0) == 0) {
                    return;
                }

                submul(basis[i ^ 1][0], basis[i ^ 0][0], a);
                submul(basis[i ^ 1][1], basis[i ^ 0][1], a);
            }
        }
    }

    void lattice_vec(
            mpz_class &r0, mpz_class &r1,
            mpz_class const &v0, mpz_class const &v1) {
        mul(r0, v0, basis[0][0]);
        addmul(r0, v1, basis[1][0]);

        mul(r1, v0, basis[0][1]);
        addmul(r1, v1, basis[1][1]);
    }

    void approx_cvp_inverse(
            mpz_class &r0, mpz_class &r1,
            mpz_class const &v0, mpz_class const &v1) {
        mul(r0, v1, basis[1][0]);
        submul(r0, v0, basis[1][1]);

        mul(r1, v0, basis[0][1]);
        submul(r1, v1, basis[0][0]);

        add(r0, half_s);
        add(r1, half_s);

        fdiv_q(r0, s);
        fdiv_q(r1, s);
    }

};

int main(int argc, char *argv[]) {
    cout << "NECA - Not Even Coppersmith's Attack" << endl;
    cout << "ROCA weak RSA key attack by Jannis Harder (me@jix.one)" << endl;
    cout << endl;
    cout << " *** Currently only 512-bit keys are supported ***" << endl;
    cout << endl;
#ifdef _OPENMP
    cout << " *** OpenMP support enabled ***" << endl;
    cout << endl;
#endif

    if (argc < 2) {
        cout << "Usage: neca <N>" << endl;
        return 0;
    }

    mpz_class n;

    if (set_str(n, argv[1], 0) != 0) {
        cout << "Could not parse RSA modulus, " <<
            "use 0x prefix for hexadecimal." << endl;
        return 1;
    }

    using namespace std;
    using std::chrono::seconds;

    std::vector<unsigned long> primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 61, 67, 71,
        73, 79, 89, 97, 101, 103, 109, 113, 127, 131, 137, 151, 157, 163
    };

    cout << "N = " << n << endl;

    neca neca_inst(n, mpz_class(65537), primes);

    if (neca_inst.dlog_n < 0) {
        cout << "Given key does not seem to be weak." << endl;
        return 1;
    } else {
        cout << "Factoring...\n" << endl;
    }

    unsigned batch_size = 10000;

    mpz_class i_start = neca_inst.dlog_n / 2;
    mpz_class i_range = neca_inst.order / 2;
    i_range /= batch_size;
    unsigned iterations = i_range.get_ui() + 1;


    unsigned completed_iterations = 0;
    bool solution_found = false;

    auto start = std::chrono::steady_clock::now();

    auto next_msg = start + seconds(1);

    #pragma omp parallel firstprivate(neca_inst)
    {
        #pragma omp for
        for (unsigned i = 0; i < iterations; ++i) {
            bool local_cancel;
            #pragma omp atomic read
            local_cancel = solution_found;

            if (local_cancel) {
                continue;
            }

            unsigned completed;
            #pragma omp atomic read
            completed = completed_iterations;

            #pragma omp critical(user_io)
            {
                auto now = std::chrono::steady_clock::now();
                if (now > next_msg) {
                    next_msg += std::chrono::seconds(1);
                    auto elapsed = now - start;


                    float progress = completed;
                    progress /= iterations;

                    auto duration = elapsed / progress;

                    auto left = duration - elapsed;

                    int n;
                    cout << " [";

                    for (int n = 1; n < 25; ++n) {
                        cout << (n < progress * 25 ? '=' : ' ');
                    }

                    cout << "] " <<
                        std::fixed << std::setw(5) << std::setprecision(2) <<
                            100 * progress << "% " <<
                        "elapsed: " << elapsed / seconds(1) << "s " <<
                        "left: " << left / seconds(1) << "s " <<
                        "total: " << duration / seconds(1) << "s" <<
                        "\e[0K\r" <<
                        flush;
                }
            }

            neca_inst.jump_guess(i_start + i * batch_size);
            for (unsigned j = 0; j < batch_size; ++j) {
                if (j > 0) {
                    neca_inst.next_guess();
                }
                if (neca_inst.search_factors()) {
                    #pragma omp critical(user_io)
                    {
                        cout << "\n\nFactorization found:" << endl;
                        cout << "N = " <<
                            neca_inst.p << " * " << neca_inst.q << endl;
                    }

                    #pragma omp atomic write
                    solution_found = true;
                    break;
                }
            }

            #pragma omp atomic
            completed_iterations += 1;
        }
    }

    if (!solution_found) {
        cout << "\nNo factorization found. :(" << endl;
        return 1;
    }
}
