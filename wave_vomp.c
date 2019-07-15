#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define WAVE_SAMPLE_NO 200
#define PI 3.14159265359

double Uniform( void );
double rand_normal( double mu, double sigma );

int main()
{
    srand((unsigned int)time(NULL));
    double wave[WAVE_SAMPLE_NO];
    int flg1 = 0;
    int flg2 = 0;
    
    // dim 2 coef calc
    int sampleno = 31;
    double m = sampleno;
    double m2 = m*m;
    double id = (1-m)/2;
    double c[100];
    for(int i=0; i<m; i++){
        double id2 = id*id;
        c[i] = (12 * m * id2 - m * (m2 - 1)) * 30 / (m2 * (m2 - 1) * (m2 - 4));
        id++;
    }

    // (1)pos_diff (2)height_diff (3)nois を回してmean(騙され量)とsigm(計測再現性)を確認する。
    // for(int k=0; k<10000; k++){ // 統計算出用loop
    // for(double pos_diff=0; pos_diff<=50; pos_diff+=2){
        // wave generate
        double peak = 0;
        double pos_diff = 20;
        double peak_hegiht = 0.5;
        for(int i=0; i<WAVE_SAMPLE_NO; i++){
            double x = (double)i;
            double rnd = rand_normal(0, 5);
            rnd = 0; // 騙され量理論値算出用
            wave[i] = 3000*exp(-pow(((x-WAVE_SAMPLE_NO/2)/22), 2)) 
                    + 3000*peak_hegiht*exp(-pow(((x-WAVE_SAMPLE_NO/2-pos_diff)/22), 2)) + rnd + 25;
            wave[i] = (int)(wave[i] + 0.5);
            if(wave[i]>peak)peak=wave[i];
            if(flg1==0)printf("#wave: %f\n",wave[i]);
        }
        flg1 = 1;

        // そのまま重心計算(従来通り、意味のない上限80%切りも行っておく…)
        double u=0;
        double d=0;
        for(int i=0; i<WAVE_SAMPLE_NO; i++){
            if(wave[i]>peak*0.8){
                u += peak*0.8*(double)(i);
                d += peak*0.8;
            }else if(wave[i]>peak*0.2){
                u += wave[i]*(double)(i);
                d += wave[i];
            }
        }
        printf("#peak1: %f\n",u/d);

        // dim 2 exec
        int ofs = (int)(1-m)/2;
        double wave_d[WAVE_SAMPLE_NO];
        double min = 0;
        int min_ind;
        for(int i=0; i<WAVE_SAMPLE_NO; i++){
            wave_d[i] = 0;
            if(i>(m-1)/2 && i<=WAVE_SAMPLE_NO-(m-1)/2){
                for(int k=0; k<m; k++){
                    wave_d[i] += c[k] * wave[i+k+ofs];
                }
            }
            if(wave_d[i]<min){
                min = wave_d[i];
                min_ind = i;
            }
            if(flg2==0)printf("#wave_d: %f\n",wave_d[i]);
        }
        flg2 = 1;
        printf("#min: %f #min_ind: %d\n",min,min_ind);

        // 二階微分の重心計算
        u = wave_d[min_ind]*(double)min_ind;
        d = wave_d[min_ind];
        int i=1;
        while(wave_d[min_ind+i] < 0 && wave_d[min_ind-i] < 0){
            u += wave_d[min_ind+i]*(double)(min_ind+i);
            u += wave_d[min_ind-i]*(double)(min_ind-i);
            d += wave_d[min_ind+i];
            d += wave_d[min_ind-i];
            i++;
        }
        // for(int i=0; i<WAVE_SAMPLE_NO; i++){
        //     if(-wave_d[i]>1){
        //         u += -wave_d[i]*(double)(i+1);
        //         d += -wave_d[i];
        //     }
        // }
        printf("#peak2: %f\n",u/d);
    // }

}



// 以下、正規分布乱数発生用関数
double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){
    double z = sqrt( -2.0*log(Uniform()) ) * sin( 2.0*PI*Uniform() );
    return mu + sigma*z;
}
