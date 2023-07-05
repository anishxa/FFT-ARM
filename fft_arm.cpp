#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.141592653589
#define N 16 // Sample Size

complex float point4_out[4];    // Radix 4 Butterfly Output
complex float point16_out[16];  // 16 Point Output Array

//Matrix multiplication for 4x4 matrix with 4x1 Matrix
complex float matmul4(complex float A[4][4], complex float B[4])
{
    for(int i = 0; i<4; i++)
    {
        point4_out[i] = 0;
    }

    for(int i = 0; i<4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            point4_out[i] += A[i][j]*B[j];
        }
    }
}

//Function to calculate the twiddle factors
complex float wnr(int r, int k,int e)
{
    if(e == 0)
        return cexpf(-I*2*PI*r*k/N);
    else
        return cexpf(I*2*PI*r*k/N);
}

// Radix 4 Butterfly
complex float rad4(complex float X[4], int k, int e)
{
    complex float mat[4][4] = {{1, 1, 1, 1},
                               {1, -I, -1, I},
                               {1, -1, 1, -1},
                               {1, I, -1, -I}};
    matmul4(mat, X);
    for(int i = 0; i<4; i++)
    {
        point4_out[i] *= wnr(i,k,e);
    }
    return 0;
}

// 16 Point FFT
complex float point16(complex float X[16], int e)
{
    complex float temp[4];
    int n = 0;
    //Splitting the 16 point sample into 4 parts of length 4
    for(int i = 0; i<4; i++)
    {
        int k = 0;
        for(int j = 0; j < 4; j++)
        {
            temp[j] = X[n+k];
            k += N/4;
        }

        // Performing the Radix 4 FFT on the split sample
        rad4(temp, i, e);

        // Setting the value to the 16 point output array
        k = 0;
        for(int j = 0; j < 4; j++)
        {
            point16_out[n+k] = point4_out[j];
            k += N/4;
        }
        n += 1;
    }

    // Stage 2
    for(int i = 0; i<4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            temp[j] = point16_out[j + 4*i];
        }

        rad4(temp, 0, e);
        for(int j = 0; j < 4; j++)
        {
            point16_out[j+ 4*i] = point4_out[j];
        }
    }
}

int main(void)
{
    complex float X[16] = {0,
                           955,
                           1023,
                           1023,
                           1023,
                           1023,
                           835,
                           831,
                           831,
                           831,
                           831,
                           838,
                           1023,
                           1023,
                           1023,
                           897};
    int e = 0;                      // FFT Mode
    point16(X, e);                  // 16 point FFT
    complex float out[16];          // Final Output Array
    complex float window[16];
    complex float window_out[16];
    float freq[16];                 // Array to store the frequency values
    float mag[16];                  // Array to store the magnitude values
    float fs = 20;                  // Sample frequency of 20Hz

    // Bit reversing the output
    int n = 0;
    for(int i = 0; i<4; i++)
    {
        int k = 0;
        for (int j = 0; j < 4; j++)
        {
            out[j + i*4] = point16_out[n + k];
            k += N/4;
        }
        n += 1;
    }

    // Printing the FFT output
    for(int i = 0; i<16; i++)
    {
        printf("%f + i%f\n", creal(out[i]), cimag(out[i]));
    }

    //Printing the frequency values
    printf("FREQUENCY\n");
    for(int i = 0; i<16; i++)
    {
        freq[i] = (float)i*fs/16;
        printf("%f ", freq[i]);

    }
    printf("\n");

    //Printing the magnitude values
    printf("MAGNITUDE\n");
    for(int i = 0; i<16; i++)
    {
        mag[i] = cabsf(out[i]);
        printf("%f ", mag[i]);
    }
    printf("\n");

    //CONVOLUTION
    double T = 1.0 / fs;    // Sampling period
    double f = 2.0;         // Frequency of the signal of interest
    double w_hanning[N];    // Hanning window values
    double w_hamming[N];    // Hamming window values

    for (int n = 0; n < N; n++) {
        double t = n * T; // Time at sample n
        // Hanning window equation
        w_hanning[n] = 0.5 * (1 - cos(2 * PI * f * t));
        // Hamming window equation
        w_hamming[n] = 0.54 - 0.46 * cos(2 * PI * f * t);
    }

    /*
    // Print the window values
    printf("Hanning window:\n");
    for (int n = 0; n < N; n++) {
        printf("w_hanning[%d] = %f\n", n, w_hanning[n]);
    }

    printf("\nHamming window:\n");
    for (int n = 0; n < N; n++) {
        printf("w_hamming[%d] = %f\n", n, w_hamming[n]);
    }
    */


    // Triangular method
    double w_triangular[N]; // Triangular window values

    for (int n = 0; n < N; n++)
    {
        w_triangular[n] = 1.0 - fabs((n - (N - 1) / 2.0) / ((N - 1) / 2.0)); // Triangular window value
    }

    /* Print the window values
    printf("Triangular window:\n");
    for (int n = 0; n < N; n++) {
        printf("w_triangular[%d] = %f\n", n, w_triangular[n]);
    }*/


    // Blackman method
    double w_blackman[N];           // Blackman window values

    double alpha_blackman = 0.16;   // Blackman window parameter

    for (int n = 0; n < N; n++)
    {
        double a0 = (1 - alpha_blackman) / 2;
        double a1 = 0.5;
        double a2 = alpha_blackman / 2;

        w_blackman[n] = a0 - a1 * cos(2 * PI * n / (N - 1)) + a2 * cos(4 * PI * n / (N - 1)); // Blackman window value
    }

    /* Print the window values
    printf("Blackman window:\n");
    for (int n = 0; n < N; n++) {
        printf("w_blackman[%d] = %f\n", n, w_blackman[n]);
    }*/


    // Turkey method
    double w_tukey[N];          // Tukey window values

    double alpha_turkey = 0.5;  // Tukey window parameter

    for (int n = 0; n < N; n++)
    {
        if (n <= alpha_turkey * (N - 1) / 2)
        {
            w_tukey[n] = 0.5 * (1 + cos(2 * PI * (n / alpha_turkey) - PI));
        }
        else if (n < (N - 1) * (1 - alpha_turkey / 2))
        {
            w_tukey[n] = 1;
        }
        else
        {
            w_tukey[n] = 0.5 * (1 + cos(2 * PI * (n - (N - 1)) / ((N - 1) * alpha_turkey) + PI));
        }
    }

    /* Print the window values
    printf("Tukey window:\n");
    for (int n = 0; n < N; n++) {
        printf("w_tukey[%d] = %f\n", n, w_tukey[n]);
    }*/
    // FFT of the window function
    // Window = w_hanning, w_tukey, w_blackman, w_triangular, w_hamming
    for (int i = 0; i < 16; i++)
    {
        window[i] = w_tukey[i];
    }

    point16(window, e);

    // Bit reversing the output
    n = 0;
    for(int i = 0; i<4; i++)
    {
        int k = 0;
        for (int j = 0; j < 4; j++)
        {
            window_out[j + i*4] = point16_out[n + k];
            k += N/4;
        }
        n += 1;
    }

    e = 1; // Setting the mode to Inverse FFT
    complex float Y[16];

    //S etting array Y as the input for the Inverse FFT
    for(int i = 0; i<16; i++)
    {
        Y[i] = out[i]*window_out[i]/16;
    }
    point16(Y, e);

    //Bit reversing the output
    n = 0;
    for(int i = 0; i<4; i++)
    {
        int k = 0;
        for (int j = 0; j < 4; j++)
        {
            out[j + i*4] = point16_out[n + k];
            k += N/4;
        }
        n += 1;
    }

    // Printing the output of Inverse FFT
    for(int i = 0; i<16; i++)
    {
        printf("%f\n", creal(out[i]));
    }

    return 0;
}
