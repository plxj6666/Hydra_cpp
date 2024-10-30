#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <gmpxx.h>
#include <NTL/GF2X.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/GF2E.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2X.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include "cryptopp/cryptlib.h"
#include "cryptopp/shake.h"
#include "cryptopp/hex.h"
#include "cryptopp/secblock.h"
#include <cryptopp/keccak.h> 

namespace CryptoPP {
    typedef unsigned char byte;
}

using CryptoPP::byte;
using CryptoPP::SHAKE128;
using CryptoPP::SecByteBlock;
using namespace std;
using namespace NTL;
using namespace CryptoPP;
using std::vector;
using std::min;
using std::ceil;
using NTL::ZZ_p;
using NTL::ZZ;
using NTL::conv;

// 定义辅助函数
ZZ mpz_class_to_ZZ(const mpz_class &mpz_val) {
    std::string str_val = mpz_val.get_str();
    return to_ZZ(str_val.c_str());
}

class Hydra {
public:
    mpz_class p;
    ZZ_p F;
    int kappa;
    int perms;
    int d;
    int Re_1, Re_2, Ri, Rh;
    Mat<ZZ_p> Me, Mi, Mh;
    vector<vec_ZZ_p> rc_b, rc_h, rc_r;

    Hydra(mpz_class p, int t, int kappa) {
        this->p = p;
        ZZ p_zz = mpz_class_to_ZZ(p); // 使用辅助转换函数
        ZZ_p::init(p_zz);

        this->F = ZZ_p(1);
        this->kappa = kappa;
        // 初始化SHAKE128
        this->perms = num_perms(t);
        this->d = get_d(p);
        tie(Re_1, Re_2, Ri, Rh) = get_rounds(p, kappa, d);

        // 构造循环矩阵
        Me.SetDims(4, 4);
        Me[0][0] = ZZ_p(3);
        Me[0][1] = ZZ_p(2);
        Me[0][2] = ZZ_p(1);
        Me[0][3] = ZZ_p(1);
        Me[1][0] = ZZ_p(1);
        Me[1][1] = ZZ_p(3);
        Me[1][2] = ZZ_p(2);
        Me[1][3] = ZZ_p(1);
        Me[2][0] = ZZ_p(1);
        Me[2][1] = ZZ_p(1);
        Me[2][2] = ZZ_p(3);
        Me[2][3] = ZZ_p(2);
        Me[3][0] = ZZ_p(2);
        Me[3][1] = ZZ_p(1);
        Me[3][2] = ZZ_p(1);
        Me[3][3] = ZZ_p(3);

        SHAKE128 shakeMatrices;
        initShake(p, "Matrices", shakeMatrices);
        // outputShake(shakeMatrices, 32); // shake正确！继续查看
        Mi = gen_mi(shakeMatrices);
        Mh = gen_mh(shakeMatrices);
        // cout << "Mi:" << Mi << endl;
        // cout << "Mh:" << Mh << endl;

        // outputShake(shakeMatrices, 32);
        // Initialize round constants
        SHAKE128 shakeConstants;
        initShake(p, "Constants", shakeConstants);
        // cout << "输出查看第一个shake" << endl;
        // outputShake(shakeConstants, 32); // 一样

        rc_b = gen_rc(Re_1 + Re_2 + Ri, 4, shakeConstants);
        // cout << rc_b[0] << endl;
        rc_h = gen_rc(Rh, 8, shakeConstants);

        SHAKE128 shakeRolling;
        initShake(p, "Rolling", shakeRolling);
        rc_r = gen_rc(perms - 1, 8, shakeRolling);
    }

    static int get_R_star(int kappa) {
        assert(kappa >= 80 && kappa <= 256);
        vector<int> R_star = {
            19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 23,
            23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26,
            26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 29, 29, 29, 29, 29,
            30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32,
            33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36, 36,
            36, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39,
            39, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42,
            43, 43, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 46, 46, 46, 46,
            46, 46, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49,
            49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52, 52, 52, 53, 53, 53,
            53, 53, 53, 54, 54, 54, 54
        };
        return R_star[kappa - 80];
    }

    static int get_round_num_head(mpz_class p, int kappa) {
        int R_star = get_R_star(kappa);

        double x0 = kappa / 24.0 + log2(12);
        double x1 = (kappa - log2(mpz_get_d(p.get_mpz_t()))) / 22.0 + log2(11);
        int R_hat = 3 + ceil(max(x0, x1));
        int R = ceil(1.25 * ceil(max({24, R_hat, R_star + 2})));
        return R;
    }

    static int get_round_num_internal(mpz_class p, int kappa, int d) {
        double x0 = kappa / 16.0;
        double x1 = (kappa - log2(mpz_get_d(p.get_mpz_t()))) / 12.0;
        int R_hat = 4 - floor(log2(d)) + ceil(max(x0, x1));
        int R = ceil(1.125 * ceil(max(kappa / 4.0 - log2(static_cast<double>(d)) + 6, static_cast<double>(R_hat))));
        return R;
    }

    static tuple<int, int, int, int> get_rounds(mpz_class p, int kappa, int d) {
        int Re_1 = 2;
        int Re_2 = 4;
        int Ri = get_round_num_internal(p, kappa, d);
        int Rh = get_round_num_head(p, kappa);
        return make_tuple(Re_1, Re_2, Ri, Rh);
    }

    static int get_d(mpz_class p) {
        for (int d = 3; d < p; ++d) {
            if (GCD(mpz_class_to_ZZ(d), mpz_class_to_ZZ(p - 1)) == 1) {
                return d;
            }
        }
        return -1;
    }

    ZZ_p field_element_from_shake(SHAKE128 &shake) {
        size_t bitlen = mpz_sizeinbase(p.get_mpz_t(), 2);
        size_t byte = (bitlen + 7) / 8;
        size_t word = (byte + 7) / 8;
        int i = 0;
        while (true) {
            vector<uint64_t> word_buf(word, 0);
            SecByteBlock buf(byte);
            shake.TruncatedFinal(buf, buf.size()); // 错误写成了sizeof(buf)
            for (size_t i = 0; i < word; i++) {
                uint64_t value = 0;
                for (size_t j = i * 8; j < min((i + 1) * 8, byte); j++) {
                    value |= (uint64_t)buf[j] << ((j - i * 8) * 8);
                }
                word_buf[i] = value;
            }
            mpz_class res = 0;
            for (auto it = word_buf.rbegin(); it != word_buf.rend(); ++it) {
                res = (res << 64) + *it;
            }
            if (res < p) {
                return conv<ZZ_p>(mpz_class_to_ZZ(res));
            }
            i += 1;
        }
    }
    
    ZZ_p field_element_from_shake_without_0(SHAKE128 &shake) {
        // outputShake(shake,32);
        while (true) {
            ZZ_p el = field_element_from_shake(shake);
            // el出问题，说明问题还是在field_element_from_shake函数
            // cout << el << endl;
            if (el != ZZ_p(0)) {
                return el;
            }
            // 重新启动SHAKE128对象，以确保生成新元素
            shake.Restart();
        }
    }

    Mat<ZZ_p> gen_mi(SHAKE128 &shake) {
        vec_ZZ_p lambda1;
        lambda1.SetLength(4);
        lambda1[0] = ZZ_p(1);
        lambda1[1] = ZZ_p(-1);
        lambda1[2] = ZZ_p(1);
        lambda1[3] = ZZ_p(-1);

        vec_ZZ_p lambda2;
        lambda2.SetLength(4);
        lambda2[0] = ZZ_p(1);
        lambda2[1] = ZZ_p(1);
        lambda2[2] = ZZ_p(-1);
        lambda2[3] = ZZ_p(-1);
        // outputShake(shake, 32);
        // shake正确
        while (true) {
            Mat<ZZ_p> M = gen_matrix(shake, 4);
            // cout << "以上是测试gen_mi功能" << endl;
            // break; 
            // 死循环,用break跳出调试
            // cout << M << endl;
            if (check_conditions(lambda1, M, 4) && check_conditions(lambda2, M, 4)) {
                return M;
            }
        }
    }


    bool check_minpoly_condition(const Mat<ZZ_p> &M, int size) {
        int max_period = 2 * size;
        Mat<ZZ_p> M_temp = M;
        for (int i = 1; i < max_period; ++i) {
            ZZ_pX minpoly = MinPolyMod(M_temp);
            if (deg(minpoly) != size || !DetIrredTest(minpoly)) {
                return false;
            }
            M_temp = M * M_temp;
        }
        return true;
    }

    bool check_conditions(const vec_ZZ_p &lambdas, const Mat<ZZ_p> &matrix, int size) {
        if (!IsInvertible(matrix) || !check_minpoly_condition(matrix, size)) {
            return false;
        }

        ZZ_p sum_ = ZZ_p(0);
        for (int j = 0; j < size; ++j) {
            ZZ_p inner = ZZ_p(0);
            for (int l = 0; l < size; ++l) {
                inner += matrix[j][l];
            }
            sum_ += lambdas[j] * inner;
        }

        if (sum_ == ZZ_p(0)) {
            return false;
        }

        for (int j = 0; j < size; ++j) {
            sum_ = ZZ_p(0);
            for (int l = 0; l < size; ++l) {
                sum_ += lambdas[l] * matrix[l][j];
            }
            if (sum_ == ZZ_p(0)) {
                return false;
            }
        }

        return true;
    }

    Mat<ZZ_p> gen_matrix(SHAKE128 &shake, int size) {
        Mat<ZZ_p> M;
        M.SetDims(size, size);

        // Initialize all elements to 1 (similar to self.F(1) in Sage)
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                M[i][j] = ZZ_p(1);
            }
        }
        // cout << M << endl;
        // outputShake(shake, 32);
        // shake没问题
        M[0][0] = field_element_from_shake_without_0(shake);
        // M[0][0]出错
        // cout << M[0][0] << endl;
        for (int i = 1; i < size; ++i) {
            M[i][0] = field_element_from_shake_without_0(shake);
            M[i][i] = field_element_from_shake_without_0(shake);
        }
        return M;
    }

    vec_ZZ_p gen_no_zero_sum(SHAKE128 &shake, int size) {
        vec_ZZ_p v;
        v.SetLength(size);
        ZZ_p sum_ = ZZ_p(0);
        while (sum_ == ZZ_p(0)) {
            sum_ = ZZ_p(0);
            for (int i = 0; i < size; ++i) {
                v[i] = field_element_from_shake_without_0(shake);
                sum_ += v[i];
            }
        }
        return v;
    }

    Mat<ZZ_p> gen_mh(SHAKE128 &shake) {
        Mat<ZZ_p> Mh;
        Mh.SetDims(8, 8);
        ZZ_p el = field_element_from_shake(shake);
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                Mh[i][j] = el;
            }
        }
        return Mh;
    }

    vector<vec_ZZ_p> gen_rc(int rounds, int size, SHAKE128 shake) {
    // 初始化返回向量
        vector<vec_ZZ_p> rc(rounds);
        cout << "计算开始" << endl;
        for (int i = 0; i < rounds; ++i) {
            // 使用默认构造函数初始化 vec_ZZ_p
            rc[i].SetLength(size);
            for (int j = 0; j < size; ++j) {
                rc[i][j] = field_element_from_shake(shake);
            }
        }
        cout << "计算结束" << endl;
        // rc错误，则问题在上面
        // cout << "rc[0][0]:" << rc[0][0] << endl;
        return rc;
    }



    bool linear_independent(const vec_ZZ_p &v1, const vec_ZZ_p &v2) {
        if (v1.length() != v2.length() || v1.length() < 2) {
            return false;
        }
        ZZ_p factor = v1[0] * inv(v2[0]);
        for (int i = 1; i < v1.length(); ++i) {
            if (v2[i] * factor != v1[i]) {
                return true;
            }
        }
        return false;
    }
    static int num_perms(int t) {
        int t_ = t / 8;
        int t__ = t % 8;

        int perms = t_;
        if (t__ > 0) {
            perms += 1;
        }
        return perms;
    }

    vec_ZZ_p non_linear_e(vec_ZZ_p state) {
        for (int i = 0; i < state.length(); ++i) {
            state[i] = power(state[i], d);
        }
        // 第一次计算正确，均为：
        // state[1000, 512, 512, 729]
        // cout << "state:" << state << endl;
        return state;
    }

    tuple<ZZ_p, ZZ_p> get_lm_dot(const vec_ZZ_p &state) {
        assert(state.length() == 4);
        ZZ_p tmp = state[0] - state[3];
        ZZ_p dot1 = tmp - state[1] + state[2];
        ZZ_p dot2 = tmp + state[1] - state[2];
        return make_tuple(dot1, dot2);
    }

    vec_ZZ_p non_linear_i(vec_ZZ_p state) {
        ZZ_p dot1, dot2;
        tie(dot1, dot2) = get_lm_dot(state);
        dot1 *= dot1;
        ZZ_p sum_ = dot1 + dot2;
        ZZ_p prod = power(sum_, 2);

        for (int i = 0; i < state.length(); ++i) {
            state[i] += prod;
        }
        return state;
    }

    vec_ZZ_p non_linear_h(vec_ZZ_p state) {
        assert(state.length() == 8);
        ZZ_p dot = state[0] + state[1] + state[2] + state[3] - state[4] - state[5] - state[6] - state[7];
        dot *= dot;

        vec_ZZ_p out;
        out.SetLength(state.length());
        for (int i = 0; i < state.length(); ++i) {
            out[i] = state[i] + dot;
        }
        return out;
    }

    tuple<vec_ZZ_p, vec_ZZ_p> non_linear_r(vec_ZZ_p y, vec_ZZ_p z) {
        ZZ_p vy, wy, wz, vz;
        tie(vy, wy) = get_lm_dot(y);
        tie(wz, vz) = get_lm_dot(z);

        ZZ_p v = vy * vz;
        ZZ_p w = wy * wz;

        for (int i = 0; i < y.length(); ++i) {
            y[i] += v;
            z[i] += w;
        }

        return make_tuple(y, z);
    }

    vec_ZZ_p R(vec_ZZ_p state, int i) {
        assert(state.length() == 8);
        assert(rc_r.size() >= i);

        if (i == 0) {
            return state;
        }

        vec_ZZ_p y = slice(state, 0, 4);
        vec_ZZ_p z = slice(state, 4, 8);

        tie(y, z) = non_linear_r(y, z);
        vec_ZZ_p y_perm = Mi * y;
        vec_ZZ_p z_perm = Mi * z;

        state = concat_vec(y_perm, z_perm) + rc_r[i - 1];
        return state;
    }

    static vec_ZZ_p concat(const vec_ZZ_p &a, const vec_ZZ_p &b) {
        // 把a,b的元素合并
        vec_ZZ_p result;
        result.SetLength(a.length() + b.length());
        for (int i = 0; i < a.length(); ++i) {
            result[i] = a[i];
        }
        for (int i = 0; i < b.length(); ++i) {
            result[a.length() + i] = b[i];
        }
        return result;
    }

    vec_ZZ_p concat_vec(const vec_ZZ_p &a, const vec_ZZ_p &b) {
        vec_ZZ_p result;
        result.SetLength(a.length() + b.length());
        for (int i = 0; i < a.length(); ++i) {
            result[i] = a[i];
        }
        for (int i = 0; i < b.length(); ++i) {
            result[a.length() + i] = b[i];
        }
        return result;
    }

    pair<vec_ZZ_p, vec_ZZ_p> permutation_b(vec_ZZ_p state) {
        vec_ZZ_p sum_;
        sum_.SetLength(4);
        clear(sum_);

        state = Me * state;
        // cout << "state" << state << endl;  // 正确


        // 错误 rc_b计算错误，说明问题在
        // cout << "rc_b[0]:" << rc_b[0] << endl; 
        for (int i = 0; i < Re_1; ++i) {
            state = non_linear_e(state);
            state = Me * state + rc_b[i];
            sum_ += state;
        }
        // 出错，说明non_linear_e有问题
        // cout << "state:" << state << endl;
        // cout << "sum_" << sum_ << endl;
        for (int i = 0; i < Ri; ++i) {
            state = non_linear_i(state);
            state = Mi * state + rc_b[i + Re_1];
            sum_ += state;
        }
        for (int i = Re_1; i < Re_1 + Re_2; ++i) {
            state = non_linear_e(state);
            state = Me * state + rc_b[i + Ri];
            if (i < Re_1 + Re_2 - 1) {
                sum_ += state;
            }
        }
        return make_pair(state, sum_);
    }

    vec_ZZ_p permutation_h(vec_ZZ_p state, const vec_ZZ_p &K) {
        for (int i = 0; i < Rh; ++i) {
            state = non_linear_h(state);
            state = Mh * state + rc_h[i];
            state += K;
        }
        return state;
    }

    vec_ZZ_p gen_ks(int t, const vec_ZZ_p &K, const vec_ZZ_p &IV, const vec_ZZ_p &N) {
        assert(IV.length() == 3);
        assert(K.length() == 4);
        vec_ZZ_p state = concat(N, IV);
        // cout << "state" << state << endl;  state正确
        assert(state.length() == 4);
        vec_ZZ_p K_vec(K);
        // cout << "K_vec:" << K_vec << endl; K_vec没问题
        // first step
        state += K_vec;
        vec_ZZ_p z;
        // cout << "state:" << state << endl; state正确，Z是才定义的，因此permutation出问题
        tie(state, z) = permutation_b(state);
        state += K_vec;
        // second step
        vec_ZZ_p K_mat = Me * K_vec;
        vec_ZZ_p K_ = concat_vec(K_vec, K_mat);
        // cout << "K_: " << K_ << endl; K_没问题
        vec_ZZ_p keystream;
        keystream.SetLength(t);
        int perm_counter = -1;

        // cout << "state" << state << endl;
        // cout << "z:" << z << endl;
        vec_ZZ_p roll = concat_vec(state, z);
        // roll 计算的有问题
        // cout << "roll:" << roll << endl;
        vec_ZZ_p perm;
        for (int i = 0; i < t; ++i) {
            int off = i % 8;
            if (off == 0) {
                perm_counter++;
                roll = R(roll, perm_counter);
                perm = permutation_h(roll, K_);
                perm += roll; // feed forward
            }
            keystream[i] = perm[off];
        }

        return keystream;
    }

    vec_ZZ_p encrypt(const vec_ZZ_p &plains, const vec_ZZ_p &K, const vec_ZZ_p &IV, const vec_ZZ_p &N) {
        int t = plains.length();
        vec_ZZ_p keystream = gen_ks(t, K, IV, N);
        // 结果不一样，说明gen_ks方法出错
        // cout << "keystream" << keystream << endl;
        vec_ZZ_p ciphers;
        ciphers.SetLength(t);
        for (int i = 0; i < t; ++i) {
            ciphers[i] = plains[i] + keystream[i];
            cout << ciphers[i] << endl;
        }
        return ciphers;
    }

    vec_ZZ_p decrypt(const vec_ZZ_p &ciphers, const vec_ZZ_p &K, const vec_ZZ_p &IV, const vec_ZZ_p &N) {
        int t = ciphers.length();
        vec_ZZ_p keystream = gen_ks(t, K, IV, N);
        vec_ZZ_p plains;
        plains.SetLength(t);
        for (int i = 0; i < t; ++i) {
            plains[i] = ciphers[i] - keystream[i];
        }
        return plains;
    }

    ZZ_p get_f() {
        return F;
    }

    void print_mat(const Mat<ZZ_p> &mat, const string &prefix) {
        cout << prefix << " = [" << endl;
        for (int i = 0; i < mat.NumRows(); ++i) {
            cout << "    [";
            for (int j = 0; j < mat.NumCols(); ++j) {
                ZZ_p el = mat[i][j];
                if (el == ZZ_p(-1)) {
                    cout << -1;
                } else if (el == ZZ_p(-2)) {
                    cout << -2;
                } else {
                    cout << el;
                }
                if (j < mat.NumCols() - 1) {
                    cout << ", ";
                }
            }
            if (i == mat.NumRows() - 1) {
                cout << "]";
            } else {
                cout << "],";
            }
            cout << endl;
        }
        cout << "]" << endl;
    }

    void print_lm_mat(Mat<ZZ_p> &mat, const string &prefix) {
        mat[0][0] -= 1;
        for (int i = 1; i < mat.NumRows(); ++i) {
            mat[i][0] -= 1;
            mat[i][i] -= 1;
        }
        print_mat(mat, prefix);
    }

    void print_conditionals() {
        print_lm_mat(Mi, "self.Mi");
        print_lm_mat(Mh, "self.Mh");
    }

private:
    static void initShake(const mpz_class &p, const string &context, SHAKE128 &shake) {
        size_t bitlen = mpz_sizeinbase(p.get_mpz_t(), 2);
        size_t num = (bitlen + 63) / 64;

        shake.Update((const byte*)"Hydra", 5);
        shake.Update((const byte*)context.data(), context.size());

        mpz_class tmp = p;
        for (size_t i = 0; i < num; i++) {
            uint64_t prime_block = tmp.get_ui();
            byte block[8];
            for (int j = 0; j < 8; j++) {
                block[j] = (prime_block >> (j * 8)) & 0xFF;
            }
            shake.Update(block, 8);
            tmp >>= 64;
        }
    }

    static void outputShake(SHAKE128 &shake, size_t length) {
        SecByteBlock buffer(length);

        shake.TruncatedFinal(buffer, buffer.size()); 
        
        string encoded;
        HexEncoder encoder(new StringSink(encoded));
        encoder.Put(buffer, buffer.size());
        encoder.MessageEnd();
        
        cout << "SHAKE128 output: " << encoded << endl;
    }

    vec_ZZ_p slice(const vec_ZZ_p &vec, long start, long end) {
        vec_ZZ_p result;
        result.SetLength(end - start);
        for (long i = start; i < end; ++i) {
            result[i - start] = vec[i];
        }
        return result;
    }

    bool IsInvertible(const Mat<ZZ_p> &matrix) {
        return (determinant(matrix) != 0);
    }

    bool DetIrredTest(const ZZ_pX &poly) {
    // 使用 NTL 库中的 IterIrredTest 函数测试多项式的不可约性
        return IterIrredTest(poly);
    }

// Function to compute the result of a polynomial evaluated at a matrix
    void CompMod(Mat<ZZ_p> &result, const ZZ_pX &P, const Mat<ZZ_p> &M) {
        long n = M.NumRows();
        long d = deg(P);

        // Initialize result to zero matrix
        result.SetDims(n, n);
        clear(result);

        // Temporary matrix to store powers of M
        Mat<ZZ_p> M_power;
        M_power.SetDims(n, n);
        ident(M_power, n);  // M_power = I

        for (long i = 0; i <= d; i++) {
            ZZ_p coef = coeff(P, i); // Get the coefficient of the polynomial P at degree i
            if (coef != 0) {
                // Add the term coeff * M_power to the result
                Mat<ZZ_p> term;
                mul(term, M_power, coef);
                add(result, result, term);
            }

            // Update M_power to M^(i+1)
            if (i < d) {
                Mat<ZZ_p> M_temp;
                mul(M_temp, M_power, M);
                M_power = M_temp;
            }
        }
    }

    // Function to check if a matrix is a zero matrix
    bool IsZero(const Mat<ZZ_p> &M) {
        long rows = M.NumRows();
        long cols = M.NumCols();
        for (long i = 0; i < rows; i++) {
            for (long j = 0; j < cols; j++) {
                if (M[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // Function to compute the characteristic polynomial using Leverrier-Faddeev algorithm
    ZZ_pX CharacteristicPolynomial(const Mat<ZZ_p> &M) {
        long n = M.NumRows();
        Mat<ZZ_p> A = M;
        Mat<ZZ_p> B;
        ZZ_pX charPoly;
        ZZ_pX tempPoly;
        ZZ_pX coeff;
        ZZ_p trace;

        // Initialize the characteristic polynomial to 1 (constant term)
        SetCoeff(charPoly, 0, ZZ_p(1));

        for (long k = 1; k <= n; k++) {
            // B = A * M
            mul(B, A, M);
            A = B;

            // Compute the trace of A
            trace = ZZ_p(0);
            for (long i = 0; i < n; i++) {
                trace += A[i][i];
            }

            // Compute the k-th coefficient of the characteristic polynomial
            SetCoeff(coeff, 0, -trace / k);
            tempPoly = charPoly;
            mul(tempPoly, tempPoly, coeff);
            charPoly += tempPoly;
        }

        // Reverse the coefficients to match the standard form
        reverse(charPoly.rep.begin(), charPoly.rep.end());

        return charPoly;
    }

    // Function to compute the minimal polynomial of a matrix
    ZZ_pX MinPolyMod(const Mat<ZZ_p> &M) {
        // Step 1: Compute the characteristic polynomial of the matrix M
        ZZ_pX charPoly = CharacteristicPolynomial(M);

        // Step 2: Initialize minPoly to charPoly as a starting point
        ZZ_pX minPoly = charPoly;

        // Step 3: Find the minimal polynomial by checking divisors of the characteristic polynomial
        for (int d = 1; d <= deg(charPoly); d++) {
            ZZ_pX divisor;
            SetCoeff(divisor, d, 1); // Set divisor as x^d
            ZZ_pX potentialMinPoly = charPoly / divisor;

            // Check if this polynomial is indeed the minimal polynomial
            Mat<ZZ_p> M_check;
            CompMod(M_check, potentialMinPoly, M); // Evaluate polynomial at matrix M
            if (IsZero(M_check)) {
                minPoly = potentialMinPoly;
                break;
            }
        }

        return minPoly;
    }
};

int main() {
    // 初始化参数
    mpz_class p("170141183460469231731687303715884105773");
    int t = 5;
    int sec = 128;

    // 创建Hydra实例
    Hydra hydra(p, t, sec);
    ZZ_p F = hydra.get_f();

    // 初始化密钥、IV和N
    vec_ZZ_p MK, IV, N;
    MK.SetLength(4);
    IV.SetLength(3);
    N.SetLength(1);

    for (int i = 0; i < 4; ++i) {
        MK[i] = ZZ_p(0);
    }
    for (int i = 0; i < 3; ++i) {
        IV[i] = ZZ_p(1);
    }
    N[0] = ZZ_p(2);

    // 初始化输入状态
    vec_ZZ_p state_in;
    state_in.SetLength(t);
    for (int i = 0; i < t; ++i) {
        state_in[i] = ZZ_p(i);
    }
    // 打印输入查看是不是一样的
    // for (int i = 0; i < state_in.length(); ++i){
    //     cout << hex << state_in[i] << endl;
    // }
    // 加密和解密操作
    vec_ZZ_p state_out = hydra.encrypt(state_in, MK, IV, N);
    vec_ZZ_p state_out2 = hydra.decrypt(state_out, MK, IV, N);

    // 验证解密结果是否与初始输入一致
    assert(state_out2 == state_in);

    /*

    cout << "state_out2:" << endl;
    for (int i = 0; i < state_out2.length(); ++i){
        cout << hex << state_out2[i] << endl;
    } */

    // cout << "state_out:" << endl;
    // 打印输出状态的每个元素的十六进制表示
    for (int i = 0; i < state_out.length(); ++i) {
        cout << hex << state_out[i] << endl;
    }

    // hydra.print_conditionals(); // 如有需要，可以取消注释以打印条件

    return 0;
}
