#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define NULL ((void *)0)



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

    fftw_complex* data = NULL;
    data = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof(fftw_complex));
    double* bytes = NULL;
    bytes = (double*)malloc((array_size_i * array_size_j) * 2 * sizeof(double));
    double *ta = NULL;
    ta = (double*)malloc(array_size_i * array_size_j * sizeof (double));
    fftw_complex* out = NULL;
    out = (fftw_complex*)malloc(array_size_i * array_size_j * sizeof(fftw_complex));

    FILE* fp;
    fopen_s(&fp, "D:\\Programming\\C++\\Qt\\CSA\\data4.bin", "rb"); // for MinGW
//    fopen("D:\\Programming\\C++\\Qt\\CSA\\data4.bin", "rb"); // for Cygwin
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
    // s_rd
    for (i = 0; i < array_size_i; ++i) {
        for (j = 0; j < array_size_j; ++j) {
            data[array_size_i * j + i] = data[array_size_i * j + i] * exp(2 * M_PI * I * fnc
                                                                          * ta[array_size_j * i + j]);
        }
    }


/*
    fftw_plan p;
    p = fftw_plan_dft_1d(1024*320, data, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
*/



    int rank = 1; /* not 2: we are computing 1d transforms */
    int n[] = {1024}; /* 1d transforms of length 10 */
    int howmany = 320;
    int idist = 1024;
    int odist = 1024;
    int istride = 1;
    int ostride = 1; /* distance between two elements in
                                  the same column */
    int *inembed = n, *onembed = n;
    fftw_plan p = fftw_plan_many_dft(rank, n,  howmany,
                                  data,  inembed,
                                  istride,  idist,
                                  out,  onembed,
                                  ostride,  odist,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    for (i = 0; i < array_size_i; ++i)
    {
        for (j = 0; j < array_size_j; ++j)
        {
            ii = array_size_i * j + i;
            //printf("i = %d  j = %d  ", i, j);
            //printf("% .10f + %.10fi\n", creal(out[ii]), cimag(out[ii]));
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
            fprintf(fp_2, "% .10f + %.10fi\n", creal(out[ii]), cimag(out[ii]));
        }
    }
    return 0;

}
