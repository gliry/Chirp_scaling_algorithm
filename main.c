#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define NULL ((void *)0)

double zeroethOrderBessel( double x )
{
    const double eps = 0.000001;
    //  initialize the series term for m=0 and the result
    double besselValue = 0;
    double term = 1;
    double m = 0;
    //  accumulate terms as long as they are significant
    while(term  > eps * besselValue)
    {
        besselValue += term;

        //  update the term
        ++m;
        term *= (x*x) / (4*m*m);
    }
    return besselValue;
}

void buildWindow( int length, double shape, double* kaiser)
{
    //  Pre-compute the shared denominator in the Kaiser equation.
    const double oneOverDenom = 1.0 / zeroethOrderBessel( shape );

    const int N = length;
    const double oneOverN = 1.0 / N;

    for (int n = 0; n < N; ++n )
    {
        const double K = (2.0 * n * oneOverN) - 1.0;
        const double arg = sqrt(1.0 - (K * K));
        kaiser[n] = zeroethOrderBessel(shape * arg) * oneOverDenom;
    }
}

int main()
{
    // Задаем основные константы
    double R_nc = 20e3; // Наклонное расстояние от центра сцены
    int Vr = 150; // Эффективная скорость радара
    double Tr = 2.5e-6; // Ширина излученного импульса
    double Kr = 20e12; // Частотная составляющая по дальности
    double f0 = 5.3e9; // Рабочая частота радара
    int BW_dop = 80; // Доплеровская полоса пропускания
    double Fr = 60e6; // Частота дискретизации по дальности
    double Fa = 200; // Частота дискретизации по азимуту
    int Naz = 1024; // Количество линий дальности (количество строк в матрице), поставить 1024
    int Nrg = 320; // Количество точек выборке по линии дальности (количество столбцов в матрице)
    double sita_r_c = (0 * 3.14) / 180; // Угол наклона луча, в радианах
    double c = 3e8; // Скорость света

    double R0 = R_nc*cos(sita_r_c); // Ближайшее наклонное расстояние к R_nc (речь о перпендикуляре)
    double Nr = Tr*Fr; // Количество точек выборки ЛЧМ сигнала по дальности
    double BW_range = Kr*Tr; // Полоса пропускания по дальности
    double lamda = c/f0; // Длина волны
    double fnc = 2*Vr*sin(sita_r_c)/lamda; // Центральная частота Допплера, рассчитанная по (4.33)
    double La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop; // Длина антенны по азимуту, согласно (4.36)
    double beta_bw = 0.886*lamda/La_real; // Луч радара 3dB
    double La = beta_bw*R0; // Длина синтезированной апертуры
    double a_sr = Fr / BW_range; // Коэффициент передискретизации по дальности
    double a_sa = Fa / BW_dop; // Коэффициент передискретизации по азимуту

    double Mamb = round(fnc/Fa); // Допплеровское размытие

    int NFFT_r = Nrg; // Длина FFT по дальности (число отсчетов по дальности)
    int NFFT_a = Naz; // Длина FFT по азимуту (число отсчетов по азимуту)

    double R_ref = R0; // Опорная цель выбирается в центре сцены и ее ближайшее расстояние равное R_ref
    double fn_ref = fnc; // Допплеровская центральная частота опорной цели

    int i = 0; // azimuth
    int j = 0; // range
    int array_size_i = 1024;
    int array_size_j = 320;
    int ii = 0;
    int read = 0;
    double D_fn_ref_Vr = 0;

    fftw_complex* data = NULL;
    data = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof(fftw_complex));

    double* bytes = NULL;
    bytes = (double*)malloc((array_size_i * array_size_j) * 2 * sizeof(double));

    double *ta = NULL;
    ta = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    double *fa = NULL;
    fa = (double*)malloc(array_size_i * sizeof (double));

    fftw_complex* out = NULL;
    out = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof(fftw_complex));

    fftw_complex* S_RD = NULL;
    S_RD = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof(fftw_complex));

    double *D_fn_Vr = NULL;
    D_fn_Vr = (double*)malloc(array_size_i * sizeof(double));

    double *D_fn_Vr_mtx = NULL;
    D_fn_Vr_mtx = (double*)malloc(array_size_i * array_size_j * sizeof(double));

    double *K_src = NULL;
    K_src = (double*)malloc(array_size_i * sizeof (double));

    double *K_src_mtx = NULL;
    K_src_mtx = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    double *Km = NULL;
    Km = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    fftw_complex * s_sc = NULL;
    s_sc = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    double * tr = NULL;
    tr = (double*)malloc(array_size_j * sizeof (double));

    double * tr_mtx = NULL;
    tr_mtx = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    fftw_complex * s_2df_1 = NULL;
    s_2df_1 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    fftw_complex * H1 = NULL;
    H1 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    double * fr_mtx = NULL;
    fr_mtx = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    double *kaiser = NULL;
    kaiser = (double*)malloc(array_size_j * sizeof (double));

    double *W_ref = NULL;
    W_ref = (double*)malloc(array_size_i * array_size_j * sizeof (double));

    fftw_complex * H1_new = NULL;
    H1_new = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    fftw_complex * s_2df_2 = NULL;
    s_2df_2 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    fftw_complex * S_RD_2 = NULL;
    S_RD_2 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    double *R0_RCMC = NULL;
    R0_RCMC = (double*)malloc(array_size_j * sizeof (double));

    fftw_complex * Haz = NULL;
    Haz = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    fftw_complex * H2 = NULL;
    H2 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));

    fftw_complex * S_RD_3 = NULL;
    S_RD_3 = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof (fftw_complex));



    FILE* fp;
    fopen_s(&fp, "D:\\Programming\\C++\\Qt\\CSA\\data4.bin", "rb"); // for MinGW
    //fopen("D:\\Programming\\C++\\Qt\\CSA\\data4.bin", "rb"); // for Cygwin
    if (!fp)
    {
        (void)printf("Unable to open file!");
        return 1;
    }

    read = fread(bytes, 8, array_size_i * array_size_j * 2, fp);
    printf("Readed %d complex numbers\n", read / 2);
    fclose(fp);

    // s_echo
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            data[ii] = bytes[ii * 2] + bytes[ii * 2 + 1] * I;
            //printf("i = %d  j = %d  ", i, j);
            //printf("% .10f + %.10fi\n", creal(data[ii]), cimag(data[ii]));
        }
    }

    // ta
    for (j = -Nrg / 2; j < (Nrg / 2) - 1; ++j)
    {
        for (i = 0; i < Naz; ++i)
        {
            //printf("\n %d  %d  ", i, j);
            ta[i] = (double) i / Fa;
            //printf("%.4f  ", ta[i]);

        }
    }

    // fa
    for (j = 0; j < NFFT_a / 2; ++j)
    {
        fa[j] = j;
        //printf("%.4f  ", fa[j]);
    }
    i = (fa[(NFFT_a / 2) - 1] + 1) * -1;
    for (j = NFFT_a / 2; j < NFFT_a; ++j)
    {
        fa[j] = i;
        i++;
        //printf("%.4f  ", fa[j]);
    }
    for (j = 0; j < NFFT_a; ++j)
    {
        fa[j] = fnc + fa[j] * (Fa / NFFT_a);
        //printf("%.4f  ", fa[j]);
    }

    // s_rd
    for (i = 0; i < array_size_i; ++i) {
        for (j = 0; j < array_size_j; ++j) {
            data[array_size_i * j + i] = data[array_size_i * j + i] * exp(2 * M_PI * I * fnc
                                                                          * ta[array_size_j * i + j]);
        }
    }

    // S_RD
    int rank = 1; /*  we are computing 1d transforms */
    int n[] = {array_size_i}; /* 1d transforms of length 1024 */
    int howmany = array_size_j;
    int idist = array_size_i;
    int odist = array_size_i;
    int istride = 1;
    int ostride = 1; /* distance between two elements in
                                  the same column */
    int *inembed = n, *onembed = n;
    fftw_plan p = fftw_plan_many_dft(rank, n,  howmany, data, inembed, istride, idist,
                                  S_RD,  onembed, ostride,  odist, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // D_fn_Vr
    for (i = 0; i < array_size_i; ++i)
    {
       D_fn_Vr[i] = sqrt(1 - (lamda * lamda * fa[i] * fa[i] / (4 * Vr * Vr)));
       //printf("j = %d i = %d %.4f \n", j, i, D_fn_Vr[array_size_i * j + i]);
    }

    // D_fn_Vr_mtx
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
           D_fn_Vr_mtx[array_size_i * j + i] = sqrt(1 - (lamda * lamda * fa[i] * fa[i] / (4 * Vr * Vr)));
           //printf("j = %d i = %d %.4f \n", j, i, D_fn_Vr[array_size_i * j + i]);
        }
    }
    D_fn_ref_Vr = sqrt(1 - lamda * lamda * fn_ref * fn_ref / (4 * Vr * Vr));

    // K_src
    for (i = 0; i <array_size_i; ++i)
    {
        K_src[i] = (2 * Vr * Vr * f0 * f0 * f0 * D_fn_Vr[i] * D_fn_Vr[i] * D_fn_Vr[i])
                / (c * R_ref * fa[i] * fa[i]);
    }

    // K_src_mtx
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            K_src_mtx[array_size_i * j + i] = (2 * Vr * Vr * f0 * f0 * f0 * D_fn_Vr[i] * D_fn_Vr[i] * D_fn_Vr[i])
                    / (c * R_ref * fa[i] * fa[i]);
            //printf("i = %d j = %d %.4f \n", i, j, K_src_mtx[array_size_i * j + i]);
        }
    }

    // Km
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            Km[array_size_i * j + i] = Kr / (1 - Kr / (K_src_mtx[array_size_i * j + i]));
        }
    }

    //tr
    for (j = 0; j < array_size_j; ++j)
    {
        tr[j] = 2 * R0 / c + (j - array_size_j / 2) / Fr;
    }
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            tr_mtx[array_size_i * j + i] = tr[j];
        }
    }

    // s_sc
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            s_sc[ii] = cexp(I * M_PI * Km[ii] * (D_fn_ref_Vr / D_fn_Vr_mtx[ii] - 1) *
                    (tr_mtx[ii] - 2 * R_ref / c * D_fn_Vr_mtx[ii]) *
                    (tr_mtx[ii] - 2 * R_ref / c * D_fn_Vr_mtx[ii]));
        }
    }

    //S_RD_1
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            S_RD[array_size_j * i + j] = S_RD[array_size_i * j + i] * s_sc[array_size_i * j + i];
        }
    }

    // s_2df_1
    int rank_2 = 1; /*  we are computing 1d transforms */
    int n_2[] = {array_size_j}; /* 1d transforms of length ... */
    int howmany_2 = array_size_i;
    int idist_2 = array_size_j;
    int odist_2 = array_size_j;
    int istride_2 = 1;
    int ostride_2 = 1; /* distance between two elements in
                                  the same column */
    int *inembed_2 = n_2, *onembed_2 = n_2;
    fftw_plan p_2 = fftw_plan_many_dft(rank_2, n_2,  howmany_2, S_RD, inembed_2, istride_2, idist_2,
                                  out,  onembed_2, ostride_2,  odist_2, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_2);
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            s_2df_1[array_size_i * j + i] = out[array_size_j * i + j];
        }
    }

    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            fr_mtx[array_size_i * j + i] = (j - NFFT_r / 2) * (Fr / NFFT_r);
        }
    }

    //
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            H1[ii] = cexp(I * M_PI * D_fn_Vr_mtx[ii] / (D_fn_ref_Vr * Km[ii]) * fr_mtx[ii] * fr_mtx[ii]) *
            cexp(I * M_PI * 4 / c * (1 / D_fn_Vr_mtx[ii] - 1 / D_fn_ref_Vr) * R_ref * fr_mtx[ii]);
        }
    }

    // (kaiser(Nrg,3).')
    buildWindow(Nrg, 3, kaiser);

    // W_ref
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            W_ref[array_size_i * j + i] = kaiser[j];
        }
    }

    // fftshift (H1,2)  swaps halves of each row
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j / 2; ++j)
        {
            H1_new[array_size_i * j + i] = H1[array_size_i * (j + array_size_j / 2) + i];
        }
    }
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = array_size_j / 2; j < array_size_j; ++j)
        {
            H1_new[array_size_i * j + i] = H1[array_size_i * (j - array_size_j / 2) + i];
        }
    }

    // S_2df_2
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            s_2df_2[array_size_j * i + j] = s_2df_1[ii] * H1_new[ii];
        }
    }

    // S_RD_2
    int rank_3 = 1; /*  we are computing 1d transforms */
    int n_3[] = {array_size_j}; /* 1d transforms of length ... */
    int howmany_3 = array_size_i;
    int idist_3 = array_size_j;
    int odist_3 = array_size_j;
    int istride_3 = 1;
    int ostride_3 = 1; /* distance between two elements in
                                  the same column */
    int *inembed_3 = n_3, *onembed_3 = n_3;
    fftw_plan p_3 = fftw_plan_many_dft(rank_3, n_3,  howmany_3, s_2df_2, inembed_3, istride_3, idist_3,
                                  out,  onembed_3, ostride_3,  odist_3, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_3);
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            S_RD_2[array_size_i * j + i] = out[array_size_j * i + j] / array_size_j;
        }
    }

    // R0_RCMC
    for (j = 0; j < array_size_j; ++j)
    {
        R0_RCMC[j] = (c / 2) * tr[j];
    }

    // Haz
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            Haz[ii] = cexp(I * 4 * M_PI * D_fn_Vr[i] * R0_RCMC[j] * f0 / c);
        }
    }

    // H2
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            H2[ii] = cexp(I * 4 * M_PI * Km[ii] / c / c * (1 - D_fn_Vr_mtx[ii] / D_fn_ref_Vr) *
             ((1 / D_fn_Vr[i]) * R0_RCMC[j] - R_ref / D_fn_Vr_mtx[ii]) *
                    ((1 / D_fn_Vr[i]) * R0_RCMC[j] - R_ref / D_fn_Vr_mtx[ii]));
        }
    }

    // S_RD_3
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            S_RD_3[ii] = S_RD_2[ii] * Haz[ii] * H2[ii];
        }
    }


    FILE* fp_2;
    fopen_s(&fp_2, "D:\\Programming\\C++\\Qt\\CSA\\data.txt", "w");
    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            fprintf(fp_2, "i = %d  j = %d  ", i, j);
            fprintf(fp_2, "% .20f + %.20fi\n", creal(S_RD_3[ii]), cimag(S_RD_3[ii]));

            //fprintf(fp_2, "i = %d j = %d %.16f \n", i, j, R0_RCMC[j]);
        }
    }
    return 0;

}





