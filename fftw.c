// FFT
// arg
//  -period X : interval default : 1
//  -interval X : interval default : 1
//  -LPF X: low path filter output filting data
//  -HPF X: high path filter ouput filting data
//  -hann    : hann window
//  -hamming : hamming window

// compile memo
// gcc fftw.c -lfftw3
// ftp://ftp.fftw.org/pub/fftw/fftw3win32mingw2.zip ⇒ a, dll, hを置けばそれで大丈夫なはず

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h> // complex.h は fftw3.h より先に include する
#include <fftw3.h>

#define MAX_LINE 200000 // 入力最大行数
#define PI 3.14159265358979323846264338327950288
#define HANN 1
#define HAMMING 2

int main(int argc, char *argv[])
{
    int i, N;
    int flg_wind = 0;
    double *x, k;
    double period = 1, lpf_ind = -1, hpf_ind = -1;

    fftw_complex *in, *out;
    fftw_plan p;

    // ****************
    // option setting
    // ****************
    if(argc > 1){
        for(i=0; i<argc; i++){
            // window set
            if(strcmp(argv[i], "-hann") == 0)
                flg_wind = HANN;
            if(strcmp(argv[i], "-hamming") == 0)
                flg_wind = HAMMING;
        }
    }
    if(argc > 2){
        for(i=0; i<argc; i++){
            // filter set
            if(strcmp(argv[i], "-lpf") == 0)
                lpf_ind = atof(argv[i+1]);
            if(strcmp(argv[i], "-hpf") == 0)
                hpf_ind = atof(argv[i+1]);
            
            // interval set
            if(strcmp(argv[i], "-period") == 0 ||
               strcmp(argv[i], "-interval") == 0)
                period = atof(argv[i+1]);

            // wind set 2
            if(strcmp(argv[i], "-wind") == 0){
                if(strcmp(argv[i+1], "hann") == 0)
                    flg_wind = HANN;
                if(strcmp(argv[i+1], "hamming") == 0)
                    flg_wind = HAMMING;
            }
        }
    }

    // **************
    // stdin
    // **************
    N = 0;
    x = (double*)malloc(sizeof(double)*MAX_LINE);
    while(fscanf(stdin, "%lf", &x[N]) != EOF){
        if(++N > MAX_LINE){
            printf("data over\n");
            free(x);
            return;
        }
    }

    k = 2 / (double)N;

    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+200));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+200));

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // **************
    //  data input 
    // **************
    for(i=0; i<N; i++){
        in[i] = x[i] + 0.0*I;
        double omg = 2*PI*(double)i/(double)N;
        if(flg_wind == HANN)
            in[i] *= 0.5 - 0.5*cos(omg);
        else if(flg_wind == HAMMING)
            in[i] *= 0.54 - 0.46*cos(omg);
    }

    // **************
    // fft exec
    // **************
    fftw_execute(p);

    // output(x axis frequency)
    if(lpf_ind < 0 && hpf_ind < 0){
        for(i=0; i<=N/2; i++){
            double re = creal(out[i]*k);
            double im = cimag(out[i]*k);
            double hz = (double)(i)/(double)N/(double)period;
            printf("%lf %g %g %g %g\n",hz, sqrt(re*re + im*im), atan2(im, re)/PI, re, im);
        }
        goto L_END;
    }

    // **************
    // filter excute
    // **************
    for(i=0; i<N; i++){
        double hz = (double)(i)/(double)N/(double)period;
        if( (lpf_ind > 0 && !(hz > lpf_ind && hz < 1/period - lpf_ind) ) || // lpf
            (hpf_ind > 0 && (hz >= hpf_ind && hz <= 1/period - hpf_ind)) ){ // hpf
                in[i] = out[i];
        }else{
            in[i] = 0.0 + 0.0*I;
        }
    }

    // **************
    // ifft exec
    // **************
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p);

    // **************
    // backfoward
    // **************
    for(i=0; i<N; i++){
        double re = creal(out[i])*k/2;
        double im = cimag(out[i])*k/2;
        double tmp = 1;
        double omg = 2*PI*(double)i/(double)N;
        if(flg_wind == HANN){
            tmp = 0.5 - 0.5*cos(omg);
        }else if(flg_wind == HAMMING){
            tmp = 0.54 - 0.46*cos(omg);
        }
        if(tmp > DBL_MIN * 100)re /= tmp;
        else re = 0;
        printf("%g\n",re);
    }

 L_END:
    // **************
    // memory release
    // **************
    if(p)fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    free(x);
    x = NULL;
    return 0;
}
