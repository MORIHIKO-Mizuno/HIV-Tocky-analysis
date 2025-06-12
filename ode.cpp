#include <R.h>
#include <math.h>
#include <vector>
#include <cmath>
#include <cstring>

extern "C"
{
    static double parms[28];
#define JQ1_w_infected parms[0]

#define JQ1_u2_infected parms[1]
#define JQ1_u3_infected parms[2]
#define JQ1_u4_infected parms[3]
#define JQ1_u5_infected parms[4]
#define JQ1_u6_infected parms[5]

#define JQ1_g_infected parms[6]

#define notreat_w_infected parms[7]

#define notreat_u2_infected parms[8]
#define notreat_u3_infected parms[9]
#define notreat_u4_infected parms[10]
#define notreat_u5_infected parms[11]
#define notreat_u6_infected parms[12]

#define notreat_g_infected parms[13]

#define stim_w_infected parms[14]

#define stim_u2_infected parms[15]
#define stim_u3_infected parms[16]
#define stim_u4_infected parms[17]
#define stim_u5_infected parms[18]
#define stim_u6_infected parms[19]

#define stim_g_infected parms[20]

#define k parms[21]
#define abp parms[22]
#define arn parms[23]
#define fb parms[24]
#define fp parms[25]
#define fr parms[26]

#define segment parms[27]

    void initparms(void (*odeparms)(int *, double *))
    {
        int N = 28;
        odeparms(&N, parms);
    }

#define xib_JQ1_infected var[0]
#define xip_JQ1_infected var[1]
#define xlr_JQ1_infected var[2]
#define xln_JQ1_infected var[3]

#define xib_notreat_infected var[4]
#define xip_notreat_infected var[5]
#define xlr_notreat_infected var[6]
#define xln_notreat_infected var[7]

#define xib_stim_infected var[8]
#define xip_stim_infected var[9]
#define xlr_stim_infected var[10]
#define xln_stim_infected var[11]

#define dxibdt_JQ1_infected vardot[0]
#define dxipdt_JQ1_infected vardot[1]
#define dxlrdt_JQ1_infected vardot[2]
#define dxlndt_JQ1_infected vardot[3]

#define dxibdt_notreat_infected vardot[4]
#define dxipdt_notreat_infected vardot[5]
#define dxlrdt_notreat_infected vardot[6]
#define dxlndt_notreat_infected vardot[7]

#define dxibdt_stim_infected vardot[8]
#define dxipdt_stim_infected vardot[9]
#define dxlrdt_stim_infected vardot[10]
#define dxlndt_stim_infected vardot[11]

#define Fb_JQ1_infected varout[0]
#define Fp_JQ1_infected varout[1]
#define Fr_JQ1_infected varout[2]

#define Fb_notreat_infected varout[3]
#define Fp_notreat_infected varout[4]
#define Fr_notreat_infected varout[5]

#define Fb_stim_infected varout[6]
#define Fp_stim_infected varout[7]
#define Fr_stim_infected varout[8]

    double delta = 0;
    // double g = 9.77061871;
    // double k = 46778.06709;
    double d = 0.377;
    double T_JQ1_uninfected, T_notreat_uninfected, T_stim_uninfected;
    double T_JQ1_infected, T_notreat_infected, T_stim_infected;

    double JQ1_uu, notreat_uu, stim_uu;
    double uuinfected = 0;
    void derivs(int *neq, double *t, double *var, double *vardot, double *varout, int *ip)
    {
        if (ip[0] < 1)
        {
            error("nout should be at least 1");
        }

        if (segment == 2)
        {
            JQ1_uu = JQ1_u2_infected;
            notreat_uu = notreat_u2_infected;
            stim_uu = stim_u2_infected;
        }
        else if (segment == 3)
        {
            JQ1_uu = JQ1_u3_infected;
            notreat_uu = notreat_u3_infected;
            stim_uu = stim_u3_infected;
        }
        else if (segment == 4)
        {
            JQ1_uu = JQ1_u4_infected;
            notreat_uu = notreat_u4_infected;
            stim_uu = stim_u4_infected;
        }
        else if (segment == 5)
        {
            JQ1_uu = JQ1_u5_infected;
            notreat_uu = notreat_u5_infected;
            stim_uu = stim_u5_infected;
        }
        else if (segment == 6)
        {
            JQ1_uu = JQ1_u6_infected;
            notreat_uu = notreat_u6_infected;
            stim_uu = stim_u6_infected;
        }
        else
        {
            fprintf(stderr, "Unknown segment: %f\n", segment);
            exit(EXIT_FAILURE);
        }

        // JQ1 infected
        T_JQ1_infected = xib_JQ1_infected + xip_JQ1_infected + xlr_JQ1_infected + xln_JQ1_infected;
        Fb_JQ1_infected = fb * xib_JQ1_infected;
        Fp_JQ1_infected = fp * xip_JQ1_infected;
        Fr_JQ1_infected = fr * xlr_JQ1_infected;
        dxlndt_JQ1_infected = arn * xlr_JQ1_infected - JQ1_uu * xln_JQ1_infected + JQ1_g_infected * xln_JQ1_infected * (1 - (T_JQ1_infected / k));
        dxibdt_JQ1_infected = JQ1_uu * xln_JQ1_infected - abp * xib_JQ1_infected - d * xib_JQ1_infected;
        dxipdt_JQ1_infected = abp * xib_JQ1_infected - JQ1_w_infected * xip_JQ1_infected - d * xip_JQ1_infected;
        dxlrdt_JQ1_infected = JQ1_w_infected * xip_JQ1_infected - arn * xlr_JQ1_infected;

        // notreat infected
        T_notreat_infected = xib_notreat_infected + xip_notreat_infected + xlr_notreat_infected + xln_notreat_infected;
        Fb_notreat_infected = fb * xib_notreat_infected;
        Fp_notreat_infected = fp * xip_notreat_infected;
        Fr_notreat_infected = fr * xlr_notreat_infected;
        dxlndt_notreat_infected = arn * xlr_notreat_infected - notreat_uu * xln_notreat_infected + notreat_g_infected * xln_notreat_infected * (1 - (T_notreat_infected / k));
        dxibdt_notreat_infected = notreat_uu * xln_notreat_infected - abp * xib_notreat_infected - d * xib_notreat_infected;
        dxipdt_notreat_infected = abp * xib_notreat_infected - notreat_w_infected * xip_notreat_infected - d * xip_notreat_infected;
        dxlrdt_notreat_infected = notreat_w_infected * xip_notreat_infected - arn * xlr_notreat_infected;

        // stim infected
        T_stim_infected = xib_stim_infected + xip_stim_infected + xlr_stim_infected + xln_stim_infected;
        Fb_stim_infected = fb * xib_stim_infected;
        Fp_stim_infected = fp * xip_stim_infected;
        Fr_stim_infected = fr * xlr_stim_infected;
        dxlndt_stim_infected = arn * xlr_stim_infected - stim_uu * xln_stim_infected + stim_g_infected * xln_stim_infected * (1 - (T_stim_infected / k));
        dxibdt_stim_infected = stim_uu * xln_stim_infected - abp * xib_stim_infected - d * xib_stim_infected;
        dxipdt_stim_infected = abp * xib_stim_infected - stim_w_infected * xip_stim_infected - d * xip_stim_infected;
        dxlrdt_stim_infected = stim_w_infected * xip_stim_infected - arn * xlr_stim_infected;
    }
}
