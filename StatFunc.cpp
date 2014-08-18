/*
 * Implementations of the statistical functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "StatFunc.h"

////////// P-value Calculatiion Functions Start ////////////////

double StatFunc::t_prob(double df, double t_value, bool two_tail) {
    double p = StatFunc::betai(df * 0.5, 0.5, df / (df + (t_value * t_value)));

    if (two_tail) return p;

    return p * 0.5;
}

double StatFunc::F_prob(double df_1, double df_2, double F_value) {
    return StatFunc::betai(df_2 * 0.5, df_1 * 0.5, df_2 / (df_2 + df_1 * F_value));
}

double StatFunc::betai(double a, double b, double x) {
    double bt;

    if (x < 0.0 || x > 1.0) throw ("Bad x in routine betai!");

    if (x == 0.0 || x == 1.0) bt = 0.0;
    else bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));

    if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
    else return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

double StatFunc::gammln(double xx) {
    int j;
    double x, y, tmp, ser;
    static const double cof[6] = {76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
        -0.5395239384953e-5};

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

double StatFunc::betacf(double a, double b, double x) {
    const int MAXIT = 300;
    const double Eps = 1.0e-08;
    const double FPMIN = (std::numeric_limits<double>::min)() / numeric_limits<double>::epsilon();
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    h = d;
    for (m = 1; m <= MAXIT; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2)*(a + m2));
        d = 1.0 + aa*d;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d*c;
        aa = -(a + m)*(qab + m) * x / ((a + m2)*(qap + m2));
        d = 1.0 + aa*d;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if (CommFunc::Abs(del - 1.0) <= Eps) break;
    }
    return h;
}

double StatFunc::chi_prob(double df, double chi_sqr_val) {
    return 1 - StatFunc::gammp(df * 0.5, chi_sqr_val * 0.5);
}

double StatFunc::gammp(const double a, const double x) {
    double gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0) throw ("Invalid arguments in routine gammp");

    if (x < a + 1.0) {
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        gcf(gammcf, a, x, gln);
        return 1.0 - gammcf;
    }
}

void StatFunc::gser(double &gamser, const double a, const double x, double &gln) {
    const int ITMAX = 500;
    const double Eps = 1.0e-08;
    int n;
    double sum, del, ap;

    gln = gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) throw ("x less than 0 in routine gser");
        gamser = 0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 0; n < ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (CommFunc::Abs(del) < CommFunc::Abs(sum) * Eps) {
                gamser = sum * exp(-x + a * log(x) - gln);
                return;
            }
        }
        throw ("a too large, ITMAX too small in routine gser");
        return;
    }
}

void StatFunc::gcf(double &gammcf, const double a, const double x, double &gln) {
    const int ITMAX = 100;
    const double EPS = numeric_limits<double>::epsilon();
    const double FPMIN = (std::numeric_limits<double>::min)() / EPS;
    int i;
    double an, b, c, d, del, h;

    gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if (CommFunc::Abs(del - 1.0) <= EPS) break;
    }
    if (i > ITMAX) throw ("a too large, ITMAX too small in gcf");
    gammcf = exp(-x + a * log(x) - gln) * h;
}
////////// P-value Calculatiion Functions End ////////////////

///////// Random Number Generation Functions Start ////////

void StatFunc::gasdev_seq(int &idum, vector<double> &vec, int size, double means, double var) {
    vec.clear();
    vec.resize(size);

    int i = 0;
    for (i = 0; i < size; i++) vec[i] = means + sqrt(var) * StatFunc::gasdev(idum);
}

void StatFunc::gasdev_seq(int &idum, vector<double> &vec, int size, double var) {
    if (size < 2) throw ("Invalid size! StatFunc::gasdev_seq");
    if (CommFunc::FloatEqual(var, 0.0)) {
        vec.clear();
        vec.resize(size);
        return;
    }

    int i = 0;
    double ave = 0.0, s_square = 0.0;
    double n = (double) size;

    while (s_square < 0.01 * var) {
        StatFunc::gasdev_seq(idum, vec, size, 0.0, var);
        for (i = 0; i < size; i++) ave += vec[i];
        ave /= n;
        for (i = 0; i < size; i++) vec[i] -= ave;
        for (i = 0, s_square = 0.0; i < size; i++) s_square += vec[i] * vec[i];
        s_square = s_square / (n - 1.0);
    }
    double c = sqrt(var / s_square);
    for (i = 0; i < size; i++) vec[i] *= c;
}

double StatFunc::gasdev(int &idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (idum < 0) iset = 0;
    if (iset == 0) {
        do {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
}

double StatFunc::UniformDev(double a, double b, int &idum) {
    if (a >= b) throw ("b must larger than a! StatFunc::UniformDev");
    if (idum > 0) idum *= -1;
    return a + (b - a) * ran1(idum);
}

double StatFunc::ran1(int &idum) {
    const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32;
    const int NDIV = (1 + (IM - 1) / NTAB);
    const double EPS = 3.0e-16, AM = 1.0 / IM, RNMX = (1.0 - EPS);
    static int iy = 0;
    static vector<int> iv(NTAB);
    int j, k;
    double temp;

    if (idum <= 0 || !iy) {
        if (-idum < 1) idum = 1;
        else idum = -idum;
        for (j = NTAB + 7; j >= 0; j--) {
            k = idum / IQ;
            idum = IA * (idum - k * IQ) - IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }
    k = idum / IQ;
    idum = IA * (idum - k * IQ) - IR*k;
    if (idum < 0) idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

double StatFunc::chidev(int &idum, const double df) {
    if (df > 2.0) return 2.0 * cheng_gamdev(idum, df * 0.5);
    else throw ("Invalid degree of freedom! StatFunc::chidev");
}

double StatFunc::cheng_gamdev(int &idum, const double alpha) {
    double u1 = 0.0, u2 = 0.0, nu = 0.0, d_buf = 0.0;
    double beta = (alpha - 1.0);
    while (1) {
        u1 = StatFunc::UniformDev(0.0, 1.0, idum);
        u2 = StatFunc::UniformDev(0.0, 1.0, idum);
        ;
        nu = (alpha - 1.0 / (6.0 * alpha)) * u1;
        nu /= (beta * u2);
        d_buf = 2.0 * (u2 - 1.0) / beta + nu + 1.0 / nu;
        if (d_buf < 2.0) return beta * nu;
        else {
            d_buf = 2.0 * log(u2) / beta - log(nu) + nu;
            if (d_buf < 1.0) return beta * nu;
        }
    }
}

int StatFunc::RandAbs(int a, int b, int &seed) { //a,bΪ������Ҷ˵㣬seedΪ����
    int rand;
    int stk = b - a + 1;
    int stl = 2;
    while (stl < stk) stl = stl + stl;
    int modul = 4 * stl;
    stk = seed;
    int sti = 1;
    while (sti <= 1) {
        stk = 5 * stk;
        stk = stk % modul;
        stl = stk / 4 + a;
        if (stl <= b) {
            rand = stl;
            sti = sti + 1;
        }
    }
    seed = stk;
    return (rand);
}

///////// Random Number Generation Functions End ////////

double StatFunc::chi_val(double df, double prob) {
    double walk = 100.0;
    double chi_val = walk;
    double way = 0.0, preway = 0.0;
    double prob_buf = chi_prob(df, chi_val);

    if (CommFunc::Abs(prob_buf - prob) < 1.0e-08) return chi_val;

    if (prob_buf > prob) {
        preway = way = 1.0;
        chi_val += walk;
    } else {
        preway = way = -1.0;
        chi_val -= walk;
    }

    while (true) {
        prob_buf = chi_prob(df, chi_val);

        if (prob_buf > prob) way = 1.0;
        else way = -1.0;

        if (CommFunc::Abs(preway - way) > 1.0e-08) walk *= 0.5;
        chi_val += walk*way;
        preway = way;

        if (walk < 1.0e-08) break;
    }

    return chi_val;
}

double StatFunc::t_val(double df, double prob) {
    double walk = 100.0;
    double t_val = walk;
    double way = 0.0, preway = 0.0;
    double prob_buf = t_prob(df, t_val, false);

    if (CommFunc::Abs(prob_buf - prob) < 1.0e-08) return t_val;

    if (prob_buf > prob) {
        preway = way = 1.0;
        t_val += walk;
    } else {
        preway = way = -1.0;
        t_val -= walk;
    }

    while (true) {
        prob_buf = t_prob(df, t_val, false);

        if (prob_buf > prob) way = 1.0;
        else way = -1.0;

        if (CommFunc::Abs(preway - way) > 1.0e-08) walk *= 0.5;
        t_val += walk*way;
        preway = way;

        if (walk < 1.0e-08) break;
    }

    return t_val;
}

double StatFunc::F_val(double df_1, double df_2, double prob) {
    double walk = 100.0;
    double F_val = walk;
    double way = 0.0, preway = 0.0;
    double prob_buf = F_prob(df_1, df_2, F_val);

    if (CommFunc::Abs(prob_buf - prob) < 1.0e-08) return F_val;

    if (prob_buf > prob) {
        preway = way = 1.0;
        F_val += walk;
    } else {
        preway = way = -1.0;
        F_val -= walk;
    }

    while (true) {
        prob_buf = F_prob(df_1, df_2, F_val);

        if (prob_buf > prob) way = 1.0;
        else way = -1.0;

        if (CommFunc::Abs(preway - way) > 1.0e-08) walk *= 0.5;
        F_val += walk*way;
        preway = way;

        if (walk < 1.0e-08) break;
    }

    return F_val;
}

double StatFunc::ControlFDR(const vector<double> &P_Value, double alpha, bool Restrict) {
    int i = 0, Size = P_Value.size();
    if (Size <= 1) throw ("Invalid size! StatFunc::ControlFDR");
    double FDR_Threshold = 0.0;
    vector<double> P_ValueBuf(Size);
    for (i = 0; i < Size; i++) P_ValueBuf[i] = P_Value[i];
    stable_sort(P_ValueBuf.begin(), P_ValueBuf.end());
    for (i = 0; i < Size; i++) {
        double d_Buf = ((double) (i + 1)) * alpha / ((double) Size);
        if (i == Size - 1) {
            if (P_ValueBuf[i] <= d_Buf || Restrict) FDR_Threshold = d_Buf;
            else FDR_Threshold = -1.0;
            break;
        }
        if (P_ValueBuf[i] <= d_Buf && P_ValueBuf[i + 1] > d_Buf) {
            FDR_Threshold = d_Buf;
            break;
        } else if (i == Size - 1) FDR_Threshold = -1.0;
    }
    return FDR_Threshold;
}

double StatFunc::ControlFDR_Zou(const vector<double> &GenePValue, double FDR) {
    int i = 0, Size = GenePValue.size();
    double PValue = 0.0;
    vector<double> Temp(Size);
    for (i = 0; i < Size; i++) Temp[i] = GenePValue[i];
    sort(Temp.begin(), Temp.end());
    for (i = 0; i < Size; i++) {
        if (Temp[i] <= (i + 1) * FDR / Size && Temp[i + 1]>(i + 1) * FDR / Size) {
            PValue = (i + 1) * FDR / Size;
            break;
        } else if (Temp[i] <= (i + 1) * FDR / Size && i == Size - 1) {
            PValue = (i + 1) * FDR / Size;
            break;
        } else if (i == Size - 1) {
            PValue = FDR;
        }
    }
    return PValue;
}

double StatFunc::ControlFDR_Storey(vector<double> &P_Value, vector<double> &Q_Value, double CrtQ, double &FDR) {
    if (P_Value.empty()) return CrtQ;

    // Initialize Lambda
    double dBuf = 0.0;
    vector<double> Lambda;
    for (dBuf = 0.0; dBuf < 0.9001; dBuf += 0.05) Lambda.push_back(dBuf);

    // Calculate Pi0
    stable_sort(P_Value.begin(), P_Value.end());
    double Pi0 = CalcuPi0(P_Value, Lambda);

    // Calculate q-value
    int i = 0, m = P_Value.size();
    double m0 = Pi0 * (double) m;
    Q_Value.clear();
    Q_Value.resize(m);
    Q_Value[m - 1] = Pi0 * P_Value[m - 1];
    for (i = m - 2; i >= 0; i--) {
        Q_Value[i] = (m0 * P_Value[i]) / (double) (i + 1);
        if (Q_Value[i] > Q_Value[i + 1]) Q_Value[i] = Q_Value[i + 1];
    }

    // Calculate FDR-adjusted critical p-value
    bool Flag = false;
    double CrtVal = 1.0;
    for (i = 0; i < m - 1; i++) {
        if (Q_Value[i] < CrtQ && Q_Value[i + 1] > CrtQ) {
            CrtVal = 0.5 * (P_Value[i] + P_Value[i + 1]);
            Flag = true;
            break;
        }
    }
    if (!Flag) CrtVal = 0.5 * (P_Value[0] + P_Value[1]);

    // Estimate FDR
    int iBuf = 0;
    for (i = 0; i < m; i++) {
        if (P_Value[i] < CrtVal) iBuf++;
    }
    FDR = (CrtVal * m0) / (double) iBuf;
    if (FDR > 1.0) FDR = 1.0;

    return CrtVal;
}

double StatFunc::CalcuPi0(vector<double> &P_Value, vector<double> &Lambda) {
    int i = 0, j = 0, P_Size = P_Value.size(), LambdaSize = Lambda.size();
    vector<double> Pi(LambdaSize), NewPi(LambdaSize);
    for (i = 0; i < LambdaSize; i++) {
        int Count = 0;
        for (j = 0; j < P_Size; j++) {
            if (P_Value[j] > Lambda[i]) Count++;
        }
        Pi[i] = Count / (P_Size * (1.0 - Lambda[i]));
    }
    double PrePi = 0.0;
    spline(Lambda, Pi, 0.0, 0.0, NewPi);
    splint(Lambda, Pi, NewPi, Lambda[LambdaSize - 1], PrePi);
    return PrePi < 1.0 ? PrePi : 1.0;
}

void StatFunc::spline(vector<double> &x, vector<double> &y, const double yp1, const double ypn, vector<double> &y2) {
    int i, k;
    double p, qn, sig, un;

    int n = y2.size();
    vector<double> u(n - 1);
    if (yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (i = 1; i < n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2]))*(ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for (k = n - 2; k >= 0; k--)
        y2[k] = y2[k] * y2[k + 1] + u[k];
}

void StatFunc::splint(vector<double> &xa, vector<double> &ya, vector<double> &y2a, const double x, double &y) {
    int k;
    double h, b, a;

    int n = xa.size();
    int klo = 0;
    int khi = n - 1;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) throw ("Bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    y = a * ya[klo] + b * ya[khi]+((a * a * a - a) * y2a[klo]
            +(b * b * b - b) * y2a[khi])*(h * h) / 6.0;
}

double StatFunc::erf(double x) {
    double a1 = 0.070523084, a2 = 0.0422820123, a3 = 0.0092705272, a4 = 0.0001520143, a5 = 0.0002765672, a6 = 0.0000430638;
    double tmp = 1 + a1 * x + a2 * pow(x, 2) + a3 * pow(x, 3) + a4 * pow(x, 4) + a5 * pow(x, 5) + a6 * pow(x, 6);
    return (1 - pow(tmp, -16));
}

// Default: upper-tail
double StatFunc::pnorm(double x) {
    double z = 0.0;
    if(x>0) z = -1.0 * x;
    else z = x;
    double sqrt2pi = 2.50662827463;
    double t0, z1, p0;
    t0 = 1 / (1 + 0.2316419 * fabs(z));
    z1 = exp(-0.5 * z * z) / sqrt2pi;
    p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
    return x >= 0 ? p0 : 1.0 - p0;
}

double StatFunc::dnorm(double x) {
    return (0.39894228 * exp(-0.5 * pow(x, 2)));
}

double StatFunc::qnorm_sub(double x, double y) {
    return (y + 0.5 * x * pow(y, 2)+(2 * pow(x, 2) + 1) * pow(y, 3) / 6 + (6 * pow(x, 3) + 7 * x) * pow(y, 4) / 12);
}

// Default: lower-tail
double StatFunc::qnorm(double p, bool upper) {
    double x = 0;
    if (upper) p = 1 - p;
    for (int i = 0; i < 4; i++) x = x + qnorm_sub(x, (pnorm(x) - p) / dnorm(x));
    return (x);
}

double StatFunc::pchisq(double x, double df) {
    if (x < 0) return -9;

    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1; // boundary function

    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

    // Check status
    if (st != 0) return -9;

    // Return p-value
    return q;
}

double StatFunc::qchisq(double q, double df) {
    if (q < 0) return -9;
    else if (q >= 1) return 0;

    double x;
    double p = 1 - q;
    int st = 0; // error variable
    int w = 2; // function variable
    double bnd = 1; // boundary function

    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

    // Check status
    if (st != 0) return -9;

    // Return p-value
    return x;
}

//#################
// functions to calculate pchisqsum 

double StatFunc::pchisqsum(double x, VectorXd lambda) {
    double pval = psadd(x, lambda);
    if (pval > 1.0) pval = psatt(x, lambda);
    return pval;
}

double StatFunc::psadd(double x, VectorXd lambda) {
double d = lambda.maxCoeff();
if (d <= 0.0) return 2.0;
lambda = lambda.array() / d;
x = x / d;

double lmin = 0.0;
double m = lambda.minCoeff();
if (m < 0.0) lmin = 0.499995 / m;
else if (x > lambda.sum()) lmin = -0.01;
else lmin = -0.5 * (double) lambda.size() / x;
double lmax = 0.499995 / lambda.maxCoeff();

double hatzeta = Brents_Kp_min_x(lambda, x, lmin, lmax, 1e-08);
if(hatzeta > lmax + 9) return 2.0;
double sign = (hatzeta < 0.0) ? -1.0 : 1.0;
double w = sign * sqrt(2 * (hatzeta * x - K(hatzeta, lambda)));
double v = hatzeta * sqrt(Kpp(hatzeta, lambda));

// debug
//cout<<"hatzeta = "<<hatzeta<<endl;
//cout<<"w = "<<w<<endl;
//cout<<"v = "<<v<<endl;


if (fabs(hatzeta) < 1e-04) return 2.0;
else return pnorm(w + log(v / w) / w);
}

double StatFunc::psatt(double x, VectorXd lambda) {
    double sum=lambda.sum();
    if(CommFunc::FloatEqual(sum,0.0)) return 2.0;
  
    double sq_sum=lambda.dot(lambda);
    double sum_sq=sum*sum;
    double a=sq_sum/sum;
    double b=sum_sq/sq_sum;
    
    return pchisq(x/a, b);
}

double StatFunc::K(double zeta, VectorXd &lambda) {
    return -0.5*(1.0 - 2.0 * zeta * lambda.array()).log().sum();
}

double StatFunc::Kp(double zeta, VectorXd &lambda) {
    return (lambda.array() / (1.0 - 2.0 * zeta * lambda.array())).sum();
}

double StatFunc::Kpp(double zeta, VectorXd &lambda) {
    return 2.0 * (lambda.array().square() / (1.0 - 2.0 * zeta * lambda.array()).array().square()).sum();
}

double StatFunc::Kp_min_x(double zeta, VectorXd &lambda, double x) {
    return Kp(zeta, lambda) - x;
}

double StatFunc::Brents_Kp_min_x(VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol) {
    double a = lowerLimit;
    double b = upperLimit;
    double c = 0;
    double d = 1.7976931348623157E+308;

    double fa = Kp_min_x(a, lambda, x);
    double fb = Kp_min_x(b, lambda, x);

    double fc = 0;
    double s = 0;
    double fs = 0;
    
    // if f(a) f(b) >= 0 then error-exit
    if (fa * fb >= 0) {
        if (fa < fb)
            return a;
        else
            return b;
    }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (fabs(fa) < fabs(fb)) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    }

    c = a;
    fc = fa;
    bool mflag = true;
    int i = 0;

    while (!(fb == 0) && (fabs(a - b) > errorTol)) {
        if ((fa != fc) && (fb != fc))
            // Inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
        else
            // Secant Rule
            s = b - fb * (b - a) / (fb - fa);

        double tmp2 = (3 * a + b) / 4;
        if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2)))) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            if ((mflag && (fabs(b - c) < errorTol)) || (!mflag && (fabs(c - d) < errorTol))) {
                s = (a + b) / 2;
                mflag = true;
            } else
                mflag = false;
        }
        fs = Kp_min_x(s, lambda, x);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        // if |f(a)| < |f(b)| then swap (a,b) end if
        if (fabs(fa) < fabs(fb)) {
            double tmp = a;
            a = b;
            b = tmp;
            tmp = fa;
            fa = fb;
            fb = tmp;
        }
        i++;
        if (i > 1000) return upperLimit+10;
    }
    return b;
}
