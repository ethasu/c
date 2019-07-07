// FFT
// arg
//  -LPF X: low path filter output filting data
//  -HPF X: high path filter ouput filting data
//  -wind -hann    : 
//        -hamming : 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h> // complex.h は fftw3.h より先に include する
#include <fftw3.h>

#define MAX_LINE 30000
#define PI 3.14159265358979323846264338327950288
#define HANN 1
#define HAMMING 2

int main(int argc, char *argv[])
{
    int i, line, hiki, N;
    int flg_wind=0, lpf_ind=0, hpf_ind=0;
    double x[MAX_LINE], y[MAX_LINE], y2[MAX_LINE], thr, k;
    double period=1;

    fftw_complex *in, *out;
    fftw_plan p;

    // option setting
    if(argc>2){
        for(i=0; i<argc; i++){
            if(strcmp(argv[i], "-wind")==0){
                if(strcmp(argv[i+1], "-hann")==0){
                    flg_wind = HANN;
                }else if(strcmp(argv[i+1], "-hamming")==0){
                    flg_wind = HAMMING;
                }
            }
            if(strcmp(argv[i], "-lpf")==0 ||
               strcmp(argv[1], "-LPF")==0){
                lpf_ind = atoi(argv[i+1]);
            }
            if(strcmp(argv[i], "-hpf")==0 ||
               strcmp(argv[1], "-hPF")==0){
                hpf_ind = atoi(argv[i+1]);
            }
            if(strcmp(argv[i], "-period")==0){
                period = atof(argv[i+1]);
            }
        }
    }

    /*標準入力からデータの読み込み*/
    line = 0;
    while(fscanf(stdin,"%lf",&x[line]) != EOF){
        line++;
    }

    N = line + 1;
    k = 2/(double)N;

    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

    /* data input */
    for(i=0; i<line; i++){
        in[i] = x[i] + 0.0*I;
        // hann window
        if(flg_wind == 1){
            in[i] *= 0.5 - 0.5*cos(2*PI*(double)i/(double)N);
        }else if(flg_wind == 2){
            // hamming window
            in[i] *= 0.54 - 0.46*cos(2*PI*(double)i/(double)N);
        }
        // tukey window
        // in[i] *= ;
    }

    // window function

    // fft exec
    fftw_execute(p);

    // output(x axis frequency)
    if(lpf_ind==0 && hpf_ind==0){
        for(i=0; i<N/2; i++){
            out[i] *= k;
            double re = creal(out[i]);
            double im = cimag(out[i]);
            printf("%lf %g %g %g %g\n",(double)i/(double)N/(double)period, sqrt(re*re + im*im), atan2(im, re)/PI, re, im);
        }
        return 0;
    }

    /* filter excute */
    memset(in, 0+0*I, sizeof(in));
    // LPF
    if(lpf_ind != 0){
        for(i=0; i<(line+1)/2; i++){
            if(i <= lpf_ind)in[i] = out[i];
        }
        in[0] *= 0.5;
    // HPF
    }else if(hpf_ind != 0){
        for(i=0; i<(line+1)/2; i++){
            if(i >= hpf_ind)in[i] = out[i];
        }
        in[0] *= 0.5;
    }else{
        printf("defaut return %s\n",argv[1]);
        return 0;
    }

    // ifft
    if(p)fftw_destroy_plan(p);
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p);

    // backfoward
    for(i=0; i<line; i++){
        double re = creal(out[i])*k;
        if(flg_wind == HANN){
            if(0.5 - 0.5*cos(2*PI*(double)i/(double)N) > __DBL_MIN__*100){
                re /= 0.5 - 0.5*cos(2*PI*(double)i/(double)N);
            }else{
                re = 0;
            }
        }else if(flg_wind == HAMMING){
            re /= 0.54 - 0.46*cos(2*PI*(double)i/(double)N);
        }
        printf("%g\n",re);
    }

    // memory release
    if(p)fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;

}