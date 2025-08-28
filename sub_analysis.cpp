#include <R.h>
#include <math.h>
#include <vector>
#include <cmath>
#include <cstring>

extern "C"
{
    static double parms[28];
#define condition1_w_infected parms[0]

#define condition1_u2_infected parms[1]
#define condition1_u3_infected parms[2]
#define condition1_u4_infected parms[3]
#define condition1_u5_infected parms[4]
#define condition1_u6_infected parms[5]

#define condition1_g_infected parms[6]

#define condition2_w_infected parms[7]

#define condition2_u2_infected parms[8]
#define condition2_u3_infected parms[9]
#define condition2_u4_infected parms[10]
#define condition2_u5_infected parms[11]
#define condition2_u6_infected parms[12]

#define condition2_g_infected parms[13]

#define condition3_w_infected parms[14]

#define condition3_u2_infected parms[15]
#define condition3_u3_infected parms[16]
#define condition3_u4_infected parms[17]
#define condition3_u5_infected parms[18]
#define condition3_u6_infected parms[19]

#define condition3_g_infected parms[20]

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

#define xib_condition1_infected var[0]
#define xip_condition1_infected var[1]
#define xlr_condition1_infected var[2]
#define xln_condition1_infected var[3]

#define xib_condition2_infected var[4]
#define xip_condition2_infected var[5]
#define xlr_condition2_infected var[6]
#define xln_condition2_infected var[7]

#define xib_condition3_infected var[8]
#define xip_condition3_infected var[9]
#define xlr_condition3_infected var[10]
#define xln_condition3_infected var[11]

#define dxibdt_condition1_infected vardot[0]
#define dxipdt_condition1_infected vardot[1]
#define dxlrdt_condition1_infected vardot[2]
#define dxlndt_condition1_infected vardot[3]

#define dxibdt_condition2_infected vardot[4]
#define dxipdt_condition2_infected vardot[5]
#define dxlrdt_condition2_infected vardot[6]
#define dxlndt_condition2_infected vardot[7]

#define dxibdt_condition3_infected vardot[8]
#define dxipdt_condition3_infected vardot[9]
#define dxlrdt_condition3_infected vardot[10]
#define dxlndt_condition3_infected vardot[11]

#define Fb_condition1_infected varout[0]
#define Fp_condition1_infected varout[1]
#define Fr_condition1_infected varout[2]

#define Fb_condition2_infected varout[3]
#define Fp_condition2_infected varout[4]
#define Fr_condition2_infected varout[5]

#define Fb_condition3_infected varout[6]
#define Fp_condition3_infected varout[7]
#define Fr_condition3_infected varout[8]

    double delta = 0;
    // double g = 9.77061871;
    // double k = 46778.06709;
    double d = 0.377;
    double T_condition1_uninfected, T_condition2_uninfected, T_condition3_uninfected;
    double T_condition1_infected, T_condition2_infected, T_condition3_infected;

    double condition1_uu, condition2_uu, condition3_uu;
    double uuinfected = 0;
    void derivs(int *neq, double *t, double *var, double *vardot, double *varout, int *ip)
    {
        if (ip[0] < 1)
        {
            error("nout should be at least 1");
        }

        if (segment == 2)
        {
            condition1_uu = condition1_u2_infected;
            condition2_uu = condition2_u2_infected;
            condition3_uu = condition3_u2_infected;
        }
        else if (segment == 3)
        {
            condition1_uu = condition1_u3_infected;
            condition2_uu = condition2_u3_infected;
            condition3_uu = condition3_u3_infected;
        }
        else if (segment == 4)
        {
            condition1_uu = condition1_u4_infected;
            condition2_uu = condition2_u4_infected;
            condition3_uu = condition3_u4_infected;
        }
        else if (segment == 5)
        {
            condition1_uu = condition1_u5_infected;
            condition2_uu = condition2_u5_infected;
            condition3_uu = condition3_u5_infected;
        }
        else if (segment == 6)
        {
            condition1_uu = condition1_u6_infected;
            condition2_uu = condition2_u6_infected;
            condition3_uu = condition3_u6_infected;
        }
        else
        {
            fprintf(stderr, "Unknown segment: %f\n", segment);
            exit(EXIT_FAILURE);
        }

        // condition1 infected
        T_condition1_infected = xib_condition1_infected + xip_condition1_infected + xlr_condition1_infected + xln_condition1_infected;
        Fb_condition1_infected = fb * xib_condition1_infected;
        Fp_condition1_infected = fp * xip_condition1_infected;
        Fr_condition1_infected = fr * xlr_condition1_infected;
        dxlndt_condition1_infected = arn * xlr_condition1_infected - condition1_uu * xln_condition1_infected + condition1_g_infected * xln_condition1_infected * (1 - (T_condition1_infected / k));
        dxibdt_condition1_infected = condition1_uu * xln_condition1_infected - abp * xib_condition1_infected - d * xib_condition1_infected;
        dxipdt_condition1_infected = abp * xib_condition1_infected - condition1_w_infected * xip_condition1_infected - d * xip_condition1_infected;
        dxlrdt_condition1_infected = condition1_w_infected * xip_condition1_infected - arn * xlr_condition1_infected;

        // condition2 infected
        T_condition2_infected = xib_condition2_infected + xip_condition2_infected + xlr_condition2_infected + xln_condition2_infected;
        Fb_condition2_infected = fb * xib_condition2_infected;
        Fp_condition2_infected = fp * xip_condition2_infected;
        Fr_condition2_infected = fr * xlr_condition2_infected;
        dxlndt_condition2_infected = arn * xlr_condition2_infected - condition2_uu * xln_condition2_infected + condition2_g_infected * xln_condition2_infected * (1 - (T_condition2_infected / k));
        dxibdt_condition2_infected = condition2_uu * xln_condition2_infected - abp * xib_condition2_infected - d * xib_condition2_infected;
        dxipdt_condition2_infected = abp * xib_condition2_infected - condition2_w_infected * xip_condition2_infected - d * xip_condition2_infected;
        dxlrdt_condition2_infected = condition2_w_infected * xip_condition2_infected - arn * xlr_condition2_infected;

        // condition3 infected
        T_condition3_infected = xib_condition3_infected + xip_condition3_infected + xlr_condition3_infected + xln_condition3_infected;
        Fb_condition3_infected = fb * xib_condition3_infected;
        Fp_condition3_infected = fp * xip_condition3_infected;
        Fr_condition3_infected = fr * xlr_condition3_infected;
        dxlndt_condition3_infected = arn * xlr_condition3_infected - condition3_uu * xln_condition3_infected + condition3_g_infected * xln_condition3_infected * (1 - (T_condition3_infected / k));
        dxibdt_condition3_infected = condition3_uu * xln_condition3_infected - abp * xib_condition3_infected - d * xib_condition3_infected;
        dxipdt_condition3_infected = abp * xib_condition3_infected - condition3_w_infected * xip_condition3_infected - d * xip_condition3_infected;
        dxlrdt_condition3_infected = condition3_w_infected * xip_condition3_infected - arn * xlr_condition3_infected;
    }
}
