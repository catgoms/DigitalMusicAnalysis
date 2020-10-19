using System;
using System.Numerics;
using System.Threading.Tasks;
using System.Threading;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];

            Parallel.For(0, wSamp, ii =>
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            });

            //for (ii = 0; ii < wSamp; ii++)
            //{
            //    double a = 2 * pi * ii / (double)wSamp;
            //    twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            //}

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];

            Parallel.For(0, nearest, kk =>
            {
                if (kk < x.Length) {
                    compX[kk] = x[kk];
                } else {
                    compX[kk] = Complex.Zero;
                }
            });

            //for (int kk = 0; kk < nearest; kk++)
            //{
            //    if (kk < x.Length)
            //    {
            //        compX[kk] = x[kk];
            //    }
            //    else
            //    {
            //        compX[kk] = Complex.Zero;
            //    }
            //}


            int cols = 2 * nearest /wSamp;

            Parallel.For(0, wSamp / 2, jj =>
            {
                timeFreqData[jj] = new float[cols];
            });

            //for (int jj = 0; jj < wSamp / 2; jj++)
            //{
            //    timeFreqData[jj] = new float[cols];
            //}

            timeFreqData = stft(compX, wSamp);
	
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            int ii = 0;
            int jj = 0;
            int kk = 0;
            int ll = 0;
            int N = x.Length;
            float fftMax = 0;
            
            float[][] Y = new float[wSamp / 2][];

            Parallel.For(0, wSamp / 2, ll =>
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            });

            //for (ll = 0; ll < wSamp / 2; ll++)
            //{
            //    Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            //}
            
            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                Parallel.For(0, wSamp, jj =>
                {
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                });

                //for (jj = 0; jj < wSamp; jj++)
                //{
                //    temp[jj] = x[ii * (wSamp / 2) + jj];
                //}
                
                tempFFT = fftTEST(temp);
                //tempFFT = fft(temp);

                for (kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }


            }

            //Parallel.For(0, (int)(2 * Math.Floor((double)N / (double)wSamp) - 1), ii =>
            //{
            //    for (kk = 0; kk < wSamp / 2; kk++) {
            //        Y[kk][ii] /= fftMax;
            //    }
            //});

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++) {
                for (kk = 0; kk < wSamp / 2; kk++) {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }

        public static int BitReverse(int n, int bits)
        {
            int reversedN = n;
            int count = bits - 1;

            n >>= 1;
            while (n > 0) {
                reversedN = (reversedN << 1) | (n & 1);
                count--;
                n >>= 1;
            }

            return ((reversedN << count) & ((1 << bits) - 1));
        }

        private Complex[] fftTEST(Complex[] x)
        {
            int bits = (int)Math.Log(x.Length, 2);

            //for (int j = 1; j < x.Length; j++) {
            //    int swapPos = BitReverse(j, bits);
            //    if (swapPos <= j) {
            //        continue;
            //    }
            //    var temp = x[j];
            //    x[j] = x[swapPos];
            //    x[swapPos] = temp;
            //}

            Parallel.For(1, x.Length, j =>
            {
                int swapPos = BitReverse(j, bits);
                if (swapPos <= j) {
                    return;
                }
                var temp = x[j];
                x[j] = x[swapPos];
                x[swapPos] = temp;
            });

            for (int N = 2; N <= x.Length; N <<= 1) {
                for (int i = 0; i < x.Length; i += N) {
                    for (int k = 0; k < N / 2; k++) {

                        int evenIndex = i + k;
                        int oddIndex = i + k + (N / 2);
                        var even = x[evenIndex];
                        var odd = x[oddIndex];

                        double term = -2 * Math.PI * k / (double)N;
                        Complex exp = new Complex(Math.Cos(term), Math.Sin(term)) * odd;

                        x[evenIndex] = even + exp;
                        x[oddIndex] = even - exp;

                    }
                }
            }

            return x;
        }
    }
}
