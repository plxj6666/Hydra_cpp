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
#include "cryptopp/pch.h"
#include "cryptopp/secblock.h"
#include <cryptopp/config.h>
#include <cryptopp/keccak.h> 
#include <sstream>
extern "C" {
    #include "KeccakHash.h"
}

#include <iomanip> // 调试使用

namespace CryptoPP {
    typedef unsigned char byte;
    extern void KeccakF1600(word64 *state);
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

// void SHAKE::Read(byte *hash, size_t size)
// {
//     CRYPTOPP_ASSERT(hash != NULLPTR);
//     ThrowIfInvalidTruncatedSize(size);

//     m_state.BytePtr()[m_counter] ^= 0x1F;
//     m_state.BytePtr()[r()-1] ^= 0x80;

//     // FIPS 202, Algorithm 8, pp 18-19.
//     while (size > 0)
//     {
//         KeccakF1600(m_state);

//         const size_t segmentLen = STDMIN(size, (size_t)BlockSize());
//         std::memcpy(hash, m_state, segmentLen);

//         hash += segmentLen;
//         size -= segmentLen;
//     }

//     Restart();
// }


// 创建一个继承自 CryptoPP::SHAKE128 的派生类
class MySHAKE128 : public CryptoPP::SHAKE128 {
public:
    // 添加自定义的 Read() 函数声明和实现
    void Read(byte *hash, size_t size) {
        CRYPTOPP_ASSERT(hash != NULLPTR);
        ThrowIfInvalidTruncatedSize(size);

        // 确保正确的填充操作（与 Read 中一致）
        m_state.BytePtr()[m_counter] ^= 0x1F;
        m_state.BytePtr()[r()-1] ^= 0x80;

        // 循环挤出数据块
        while (size > 0)
        {
            KeccakF1600(m_state);  // 调用 Keccak 的核心函数

            const size_t segmentLen = CryptoPP::STDMIN(size, (size_t)BlockSize());
            std::memcpy(hash, m_state, segmentLen);  // 将状态中的数据复制到输出中

            hash += segmentLen;
            size -= segmentLen;
        }
    }
};


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
        // 初始化MySHAKE128
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

        Keccak_HashInstance shakeMatrices;
        if (Keccak_HashInitialize_SHAKE128(&shakeMatrices) != KECCAK_SUCCESS) {
            printf("Failed to initialize shakeMatrices\n");
        }
        initShake(p, "Matrices", shakeMatrices);
        // cout << "Matrices initialization end" << endl;
        // outputShake(shakeMatrices, 32); // shake正确！继续查看



        // cout << "start gen mi" << endl;
        Mi = gen_mi(shakeMatrices);
        // cout << "end gen mi" << endl;
        // cout << "Mi:" << Mi << endl;
        // Mi right


        Mh = gen_mh(shakeMatrices);
        // cout << "Mh:" << Mh << endl;
        // mh right

        // outputShake(shakeMatrices, 32);
        // Initialize round constants
        Keccak_HashInstance shakeConstants;
         if (Keccak_HashInitialize_SHAKE128(&shakeConstants) != KECCAK_SUCCESS) {
            printf("Failed to initialize shakeConstants\n");
        }
        initShake(p, "Constants", shakeConstants);
        // cout << "Constants initialization end" << endl;

        cout << "begin rc_b" << endl;
        rc_b = gen_rc(Re_1 + Re_2 + Ri, 4, shakeConstants);
        
        // cout << "rc_b:" << endl;
        // for (int i = 0; i < rc_b.size(); i ++) {
        //     cout << rc_b[i] << ' ';
        // }
        //rc_b right
        cout << "end rc_b" << endl; 

        cout << Rh << endl;
        rc_h = gen_rc(Rh, 8, shakeConstants);
        // cout << "rc_h:" << endl;
        // for (int i = 0; i < rc_h.size(); i ++) {
        //     cout << rc_h[i] << ' ';
        // }
        //  rc_h right
        Keccak_HashInstance shakeRolling;
         if (Keccak_HashInitialize_SHAKE128(&shakeRolling) != KECCAK_SUCCESS) {
            printf("Failed to initialize shakeRolling\n");
        }
        initShake(p, "Rolling", shakeRolling);
        // cout << "Rolling initialization end" << endl;
        // rc_r = gen_rc(perms - 1, 8, shakeRolling);
        // cout << "rc_r:" << endl;
        // for (int i = 0; i < rc_r.size(); i ++) {
        //     cout << rc_r[i] << ' ';
        // }
        // rc_r right 均为空
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

    ZZ_p field_element_from_shake(Keccak_HashInstance &shake) {
        size_t bitlen = mpz_sizeinbase(p.get_mpz_t(), 2);
        size_t byte = (bitlen + 7) / 8;
        size_t word = (byte + 7) / 8;
        int i = 0;
        while (true) {
            vector<uint64_t> word_buf(word, 0);
            SecByteBlock buf(byte);
            if (Keccak_HashSqueeze(&shake, buf, buf.size() * 8) != KECCAK_SUCCESS) {
                printf("Failed to squeeze the buf\n");
            }
            // string encoded;
            // HexEncoder encoder(new StringSink(encoded));
            // encoder.Put(buf, buf.size());
            // encoder.MessageEnd();

        
            // cout << "Test SHAKE128 output: " << encoded << endl;
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
            // cout << "res :" << res << endl;
            if (res < p) {
                // string encoded;
                // HexEncoder encoder(new StringSink(encoded));
                // encoder.Put(buf, buf.size());
                // encoder.MessageEnd();
        
                // cout << "Final SHAKE128 output: " << encoded << endl;
                return conv<ZZ_p>(mpz_class_to_ZZ(res));
            }
            i += 1;
        }
    }
    
    ZZ_p field_element_from_shake_without_0(Keccak_HashInstance &shake) {
        while (true) {
            ZZ_p el = field_element_from_shake(shake);
            if (el != ZZ_p(0)) {
                return el;
            }
        }
    }

    Mat<ZZ_p> gen_mi(Keccak_HashInstance &shake) {
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
        while (true) {
            Mat<ZZ_p> M = gen_matrix(shake, 4);
            // 问题2：由于shake的随机性，导致这个M计算不能和sage相比较，因此问题应该出在下面这个返回条件
            if (check_conditions(lambda1, M, 4) && check_conditions(lambda2, M, 4)) {
                // cout << "return successfully!" << endl;
                return M;
            }
        }
    }


    bool check_minpoly_condition(const Mat<ZZ_p> &M, int size) {
        int max_period = 2 * size;
        Mat<ZZ_p> M_temp = M;
        for (int i = 1; i < max_period; ++i) {
            ZZ_pX charPoly = CharacteristicPolynomial(M_temp);
            ZZ_pX minpoly = MinimalPolynomialFromCharPoly(charPoly);
            // cout << "Minimal Polynomial: " << minpoly << endl;
            // cout << "deg of minpoly:" << deg(minpoly) << endl;
            // cout << "isIrreducible? :" <<  DetIrredTest(minpoly) << endl;
            // both of above are wrong!
            if (deg(minpoly) != size || !DetIrredTest(minpoly)) {
                return false;
            }
            M_temp = M * M_temp;
        }
        return true;
    }

    bool check_conditions(const vec_ZZ_p &lambdas, const Mat<ZZ_p> &matrix, int size) {
        bool isInvertible = IsInvertible(matrix);
        bool is_minpoly = check_minpoly_condition(matrix, size);
        // cout << "M:" << matrix << endl;
        // cout << "isInvertible:" << isInvertible << endl;
        // cout << "isMinpoly:" << is_minpoly << endl;
        // cout << endl;
        // 假如取反
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

    Mat<ZZ_p> gen_matrix(Keccak_HashInstance &shake, int size) {
        Mat<ZZ_p> M;
        M.SetDims(size, size);

        // Initialize all elements to 1 (similar to self.F(1) in Sage)
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                M[i][j] = ZZ_p(1);
            }
        }
        M[0][0] = field_element_from_shake_without_0(shake);
        // M[0][0]不同
        // cout <<"M[0][0]:" << M[0][0] << endl;
        for (int i = 1; i < size; ++i) {
            M[i][0] = field_element_from_shake_without_0(shake);
            M[i][i] = field_element_from_shake_without_0(shake);
        }
        return M;
    }

    vec_ZZ_p gen_no_zero_sum(Keccak_HashInstance &shake, int size) {
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

    Mat<ZZ_p> gen_mh(Keccak_HashInstance &shake) {
        Mh.SetDims(8, 8);
        vec_ZZ_p lambdas;
        lambdas.SetLength(8);
        lambdas[0] = ZZ_p(1);
        lambdas[1] = ZZ_p(1);
        lambdas[2] = ZZ_p(1);
        lambdas[3] = ZZ_p(1);
        lambdas[4] = ZZ_p(-1);
        lambdas[5] = ZZ_p(-1);
        lambdas[6] = ZZ_p(-1);
        lambdas[7] = ZZ_p(-1);
        while (true) {
            Mat<ZZ_p> M = gen_matrix(shake, 8);
            if (check_conditions(lambdas, M, 8)) {
                // cout << "return successfully!" << endl;
                return M;
            }
        }
    }

    vector<vec_ZZ_p> gen_rc(int rounds, int size, Keccak_HashInstance &shake) {
    // 初始化返回向量
        vector<vec_ZZ_p> rc(rounds);
        for (int i = 0; i < rounds; ++i) {
            // 使用默认构造函数初始化 vec_ZZ_p
            rc[i].SetLength(size);
            for (int j = 0; j < size; ++j) {
                rc[i][j] = field_element_from_shake(shake);
            }
        }
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
        vec_ZZ_p ciphers;
        ciphers.SetLength(t);
        for (int i = 0; i < t; ++i) {
            ciphers[i] = plains[i] + keystream[i];
            // cout << ciphers[i] << endl;
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
    static void initShake(const mpz_class &p, const string &context, Keccak_HashInstance &shake) {
        size_t bitlen = mpz_sizeinbase(p.get_mpz_t(), 2);
        size_t num = (bitlen + 63) / 64;

        // shake.Update((const byte*)"Hydra", 5);
        if (Keccak_HashUpdate(&shake, (const byte*)"Hydra", 5 * 8) != KECCAK_SUCCESS) {
            printf("Failed to update Hydra\n");
        }

        // shake.Update((const byte*)context.data(), context.size());
        if (Keccak_HashUpdate(&shake, (const byte*)context.data(), context.size() * 8) != KECCAK_SUCCESS) {
            printf("Failed to update context\n");
        }

        mpz_class tmp = p;
        for (size_t i = 0; i < num; i++) {
            uint64_t prime_block = tmp.get_ui();
            byte block[8];
            for (int j = 0; j < 8; j++) {
                block[j] = (prime_block >> (j * 8)) & 0xFF;
            }
            // shake.Update(block, 8);
            if (Keccak_HashUpdate(&shake, block, 8 * 8) != KECCAK_SUCCESS) {
                printf("Failed to update context\n");
            }
            tmp >>= 64;
        }

            // 最终化吸收阶段
        if (Keccak_HashFinal(&shake, NULL) != KECCAK_SUCCESS) {
            printf("Failed to finalize the input data\n");
        }
    }

    // static void outputShake(Keccak_HashInstance &shake, size_t length) {
    //     SecByteBlock buffer(length);

    //     // shake.Read(buffer, buffer.size()); 
    //     if (Keccak_HashSqueeze(&shake, buffer, buffer.size() * 8) != KECCAK_SUCCESS) {
    //         printf("Failed to squeeze the buffer\n");
    //     }
        
    //     string encoded;
    //     HexEncoder encoder(new StringSink(encoded));
    //     encoder.Put(buffer, buffer.size());
    //     encoder.MessageEnd();
        
    //     cout << "MySHAKE128 output: " << encoded << endl;
    // }

    // string outputshake(SecByteBlock buf){
    //     string encoded;
    //     HexEncoder encoder(new StringSink(encoded));
    //     encoder.Put(buf, buf.size());
    //     encoder.MessageEnd();
    // }

    vec_ZZ_p slice(const vec_ZZ_p &vec, long start, long end) {
        vec_ZZ_p result;
        result.SetLength(end - start);
        for (long i = start; i < end; ++i) {
            result[i - start] = vec[i];
        }
        return result;
    }

    bool IsInvertible(const Mat<ZZ_p> &matrix) {
        return (determinant(matrix) != 0); // 返回可逆
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

    // 计算 n-1 阶子矩阵，用于行列式的递归计算
    Mat<ZZ_pX> SubMatrix(const Mat<ZZ_pX>& M, long skipRow, long skipCol) {
        long n = M.NumRows();
        Mat<ZZ_pX> subM;
        subM.SetDims(n - 1, n - 1);

        for (long i = 0, subI = 0; i < n; i++) {
            if (i == skipRow) continue;  // 跳过指定行

            for (long j = 0, subJ = 0; j < n; j++) {
                if (j == skipCol) continue;  // 跳过指定列

                subM[subI][subJ] = M[i][j];
                subJ++;
            }
            subI++;
        }

        return subM;
    }

    // 递归计算矩阵的行列式
    ZZ_pX Determinant(const Mat<ZZ_pX>& M) {
        long n = M.NumRows();
        
        if (n == 1) {
            return M[0][0];  // 对于 1x1 矩阵，行列式就是元素本身
        }

        ZZ_pX det;  // 行列式的值
        clear(det);  // 初始化为 0

        for (long j = 0; j < n; j++) {
            // 计算 (-1)^(0+j) * M[0][j] * 子矩阵的行列式
            ZZ_pX cofactor = M[0][j] * Determinant(SubMatrix(M, 0, j));

            // 根据符号交替相加减
            if (j % 2 == 0) {
                det += cofactor;  // 正符号
            } else {
                det -= cofactor;  // 负符号
            }
        }

        return det;
    }

    // 确保多项式系数严格遵循模数
    ZZ_pX ModReduce(const ZZ_pX& poly) {
        ZZ_pX reducedPoly = poly;
        for (long i = 0; i <= deg(poly); i++) {
            SetCoeff(reducedPoly, i, coeff(poly, i));  // 去掉 % ZZ_p::modulus()
        }
        return reducedPoly;
    }

    // 计算特征多项式: det(M - λI)
    ZZ_pX CharacteristicPolynomial(const Mat<ZZ_p>& M) {
        long n = M.NumRows();
        ZZ_pX lambda;  // 用来表示特征多项式中的 λ 变量
        SetX(lambda);  // 设置 lambda 为 X，即特征多项式中的未知数

        // 构造矩阵 M - λI
        Mat<ZZ_pX> A_lambda;
        A_lambda.SetDims(n, n);

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < n; j++) {
                A_lambda[i][j] = conv<ZZ_pX>(M[i][j]);  // 转换矩阵元素为多项式形式
                if (i == j) {
                    A_lambda[i][j] -= lambda;  // 对角线元素减去 λ
                }
            }
        }

        // 计算行列式，这就是特征多项式
        ZZ_pX charPoly = Determinant(A_lambda);
        return ModReduce(charPoly);  // 返回确保模数范围内的特征多项式
    }

    // 找到多项式的最小多项式，去掉重复因子
    ZZ_pX MinimalPolynomialFromCharPoly(const ZZ_pX& charPoly) {
        ZZ_pX minPoly = charPoly;
        ZZ_pX derivedPoly = diff(charPoly);  // 计算特征多项式的导数
        ZZ_pX gcdPoly;

        while ((gcdPoly = GCD(minPoly, derivedPoly)) != 1) {  // 找到最小公因式为 1 时的多项式
            minPoly /= gcdPoly;  // 消去重复因子
            derivedPoly = diff(minPoly);  // 继续求导数
        }

        return ModReduce(minPoly);  // 返回最终的最小多项式，并确保模数范围
    }
};

// 将 NTL 的 ZZ 类型转换为 GMP 的 mpz_t 类型
void ZZ_to_mpz(mpz_t& gmp_num, const ZZ& num) {
    // 使用 stringstream 将 ZZ 类型转换为字符串
    stringstream ss;
    ss << num;  // 将 ZZ 写入流中
    string str = ss.str();  // 从流中获取字符串
    // 将字符串转换为 GMP 的 mpz_t
    mpz_set_str(gmp_num, str.c_str(), 10); // 10 表示十进制
}

void printHex(const Vec<ZZ_p>& state_out) {
    for (size_t i = 0; i < state_out.length(); i++) {
        ZZ num = rep(state_out[i]);  // 获取 ZZ 类型的整数
        
        // cout << "Decimal: " << num << endl;
        
        // 将 ZZ 转换为 GMP 的 mpz_t 类型
        mpz_t gmp_num;
        mpz_init(gmp_num);
        ZZ_to_mpz(gmp_num, num);  // 将 NTL 的 ZZ 转换为 GMP 的 mpz_t

        // 使用 GMP 函数将数值转换为16进制字符串
        char* hex_str = mpz_get_str(NULL, 16, gmp_num);  // 16 表示16进制

        // 输出结果
        cout << "0x" <<  hex_str << endl;

        // 释放内存
        mpz_clear(gmp_num);
        free(hex_str);  // 释放由 GMP 动态分配的字符串
    }
}

int main() {
    // 初始化参数
    mpz_class p("170141183460469231731687303715884105773");
    int t = 5;
    int sec = 128;

    cout << "Initialization Start" << endl;
    // 创建Hydra实例
    Hydra hydra(p, t, sec);
    cout << "Initialization End" << endl;
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
    cout << state_in << endl;
    cout << MK << endl;
    cout << IV << endl;
    cout << N << endl;
    vec_ZZ_p state_out = hydra.encrypt(state_in, MK, IV, N);
    vec_ZZ_p state_out2 = hydra.decrypt(state_out, MK, IV, N);

    // 验证解密结果是否与初始输入一致
    assert(state_out2 == state_in);
    printHex(state_out);
    // hydra.print_conditionals(); // 如有需要，可以取消注释以打印条件

    return 0;
}
