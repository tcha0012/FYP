using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using DSP;
using Accord.Math;


namespace PCG_segmentation
{
    class PCG_methods
    {

        public static void csv_to_list(string path, List<int> indexes, List<double> timings, List<double> ECG, List<double> PCG, double length, int Fs)
        {
            using (var reader = new StreamReader(path))
            {
                
                int i = 0;
                bool first_row_flag = true;
                while ((!reader.EndOfStream) && (i <= Convert.ToInt32(Convert.ToDouble(Fs) * length)))
                {
                    var line = reader.ReadLine();
                    var values = line.Split(';', ',');

                    if (!first_row_flag) {
                        indexes.Add(Convert.ToInt32(values[0]));
                        timings.Add(Convert.ToDouble(values[1]));
                        ECG.Add(Convert.ToDouble(values[2]));
                        PCG.Add(Convert.ToDouble(values[3]));
                    }
                    first_row_flag = false;
                    i++;
                }
            }

        }

        public static void save_to_csv(List<double> data, string file_location)
        {
            StringBuilder to_csv = new StringBuilder();

            try
            {
                System.IO.File.Delete(file_location);
            }
            catch
            {
                return;
            }

            for (int i = 0; i < data.Count; i++)
            {
                to_csv.AppendLine(Convert.ToString(data[i]));

            }
            try
            {
                File.AppendAllText(file_location, to_csv.ToString());
            }
            catch (Exception ex)
            {
                Console.WriteLine("Data could not be written to the CSV file.");
                return;
            }
        }
        public static void get_hr_preprocessing(List<double> processed_data, List<double> audiofile, int Fs, int Fs_new)
        {

            // Matlab script reproducing:
            /*
            % dealing with padded zeros and instances of zeros -- If there is a value of 0, replace it with lowest absolute value that is not 0
            if ~isempty(min(abs(audio_data(audio_data~= 0))))
                audio_data(audio_data == 0) = min(abs(audio_data(audio_data~= 0)));
            end
            */

            // get absolute list
            save_to_csv(audiofile, @"C:\github\PCG_Segmentation\raw_data.csv");
            List<double> absolute_audio_file = new List<double>();

            for (int i = 0; i < audiofile.Count; i++)
            {
                if (audiofile[i] != 0) {
                    absolute_audio_file.Add(Math.Abs(audiofile[i]));
                }
            }
            // Finding absolute minimum value
            double minimum = absolute_audio_file.Min();
            // Replacing 0s with absolute minimum value
            for (int i = 0; i < audiofile.Count; i++)
            {
                if (audiofile[i] == 0)
                {
                    audiofile[i] = minimum;
                }
            }

            // Resampling - To simplify the process, a simplified resampling approach has been taken,
            // Assuming that the new sampling rate is a factor of the original.
            List<double> audiofile_1000 = new List<double>();
            for (int i = 0; i < audiofile.Count; i++)
            {
                if (i % (Fs/Fs_new) == 0)
                {
                    audiofile_1000.Add(audiofile[i]);
                }
            }



            ///// TESTING THE MATHNET LIBRARIES

            //(byte n, double wc1, double wc2) = MathNet.Filtering.Butterworth.Designer.BandPass(20/1000, 100/1000, 300/1000, 405/1000, 0.05, 50);

            //double lowStopbandFreq = 20 / 1000;
            //double lowPassbandFreq = 100 / 1000;
            //double highPassbandFreq = 300 / 1000;
            //double highStopbandFreq = 405 / 1000;
            //double passbandRipple = 0.05;
            //double stopbandAttenuation = 50;
            //ValueTuple<Double[], Double[]> coefficients = MathNet.Filtering.Butterworth.IirCoefficients.BandPass(lowStopbandFreq, lowPassbandFreq, highPassbandFreq, highStopbandFreq, passbandRipple, stopbandAttenuation);
            ////(byte n, double wc1, double wc2) = MathNet.Filtering.Butterworth.Designer.BandPass(lowStopbandFreq, lowPassbandFreq, highPassbandFreq, highStopbandFreq, passbandRipple, stopbandAttenuation);
            //// Performing High and low pass filter
            LowpassFilterButterworthImplementation filter_low = new LowpassFilterButterworthImplementation(400, 2, Fs_new);
            HighpassFilterButterworthImplementation filter_high = new HighpassFilterButterworthImplementation(25, 2, Fs_new);


            //float resonance = 1.414f;
            //FilterButterworth filter = new FilterButterworth(25, 1000, FilterButterworth.PassType.Highpass,resonance);

            List<double> filtered_audiofile = new List<double>();
            //double curr_data_point = 0;
            //float data_point_F = 0f;
            //for (int i = 0; i < audiofile_1000.Count(); i++)
            //{
            //    curr_data_point = audiofile_1000[i];
            //    data_point_F = (float)curr_data_point;
            //    filter.Update(data_point_F)
            //    filtered_audiofile.Add(filter.Value());
            //}

            for (int i = 0; i < audiofile_1000.Count; i++)
            {
                filtered_audiofile.Add(filter_low.compute(audiofile_1000[i]));
            }

            for (int i = 0; i < audiofile_1000.Count; i++)
            {
                filtered_audiofile[i] = filter_high.compute(filtered_audiofile[i]);
            }

            // save_to_csv(filtered_audiofile, @"C:\github\PCG_Segmentation\filtered_data.csv");
            // Viterbi_Springer.PlotData(filtered_audiofile, 1000, "PCG_amplitude", "time (s)", "PCG with 25-400 Hz bandpass filter", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PCG_bandpass.png");


            // Spike removal
            Schmidt_spike_removal.schmidt_spike_removal(filtered_audiofile, processed_data, Fs_new);

            // copying to csv so I can analyse plot.
            // save_to_csv(processed_data, @"C:\github\PCG_Segmentation\spikes_removed.csv");
        }

        static double standard_deviation(List<double> signal)
        {
            double mean_of_signal = signal.Average();
            double sum = 0;
            for (int i = 0; i<signal.Count(); i++)
            {
                sum += Math.Pow((signal[i] - mean_of_signal),2);
            }
            sum = Math.Sqrt(sum / signal.Count());
            return sum;
        }
        

        public static void normalise_signal(List<double> signal)
        {
            double mean_of_signal = signal.Average();
            double std_of_signal = standard_deviation(signal);
            for (int i = 0; i < signal.Count(); i++)
            {
                signal[i] = (signal[i] - mean_of_signal) / std_of_signal;
            }

        }

        /*
        public static void getHilbert_envelope(List<double> hilbert_envelope, List<double> audio_data, int Fs)
        {

            double[] audio_array = audio_data.ToArray();
            uint length = (uint)audio_data.Count();
            DFT dft = new DFT();

            dft.Initialize(length);
            List<Complex> audio_complex = new List<Complex>();
            for (int i = 0; i < audio_data.Count(); i++)
            {
                audio_complex.Add(new Complex(audio_data[i],0));
            }
            Complex[] audio_complex_array = audio_complex.ToArray();
            MathNet.Numerics.IntegralTransforms.Fourier.Forward(audio_complex_array);
            // Get DFT of audiodata
            //Complex[] cSpectrum = dft.Execute(audio_array);

            int n = (int)length;
            List<int> h = new List<int>();
            if ((n > 0) && (n%2 == 0)) // If even
            {
                for (int i = 0; i < n; i++)
                {
                    if (i % (n/2) == 0)
                    {
                        h.Add(1);
                    }
                    else if ((i > 1) && (i < n/2))
                    {
                        h.Add(2);
                    }
                    else
                    {
                        h.Add(0);
                    }
                }
            }
            else // if odd
            {
                for (int i = 0; i < n; i++)
                {
                    if (i == 0)
                    {
                        h.Add(1);
                    }
                    else if ((i > 1) && (i <= n / 2))
                    {
                        h.Add(2);
                    }
                    else
                    {
                        h.Add(0);
                    }
                }
            }


            // Multiply by hilbert transform

            for (int i = 0; i < h.Count(); i++)
            {
                audio_complex_array[i] = Complex.Multiply(audio_complex_array[i], (new Complex(h[i], 0)));
            }

            MathNet.Numerics.IntegralTransforms.Fourier.Inverse(audio_complex_array);
            for (int i = 0; i < n; i++)
            {
                hilbert_envelope.Add(audio_complex_array[i].Magnitude);
            }
            save_to_csv(hilbert_envelope, @"C:\github\PCG_Segmentation\hilbert_envelope.csv");

            
            normalise_signal(hilbert_envelope);
            save_to_csv(hilbert_envelope, @"C:\github\PCG_Segmentation\normalised_hilbert_envelope.csv");
            
        }
        */
        /*
        public static void fft(List<Complex> fft_input, List<Complex> fft_output, int fft_length)
        {
            for (int i =0; i < fft_length; i++)
            {
                fft_output.Add(new Complex(0, 0));
                Complex minus_i = new Complex(0, -1);
                for (int j = 0; j < fft_input.Count(); j++)
                {
                    fft_output[i] = Complex.Add(fft_output[i], Complex.Multiply(fft_input[j], Complex.Exp(Complex.Multiply(-minus_i, (2 * Math.PI * i * j) / fft_input.Count))));

                }
            }
        }
        */
        /*
        public static void get_hr_autocorrelation_filtered(List<double> signal_autocorrelation_filtered, List<double> envelope, int Fs)
        {
            List<Complex> y = new List<Complex>();
            double mean = envelope.Average();

            for(int i =0; i < envelope.Count(); i++)
            {
                y.Add(new Complex(envelope[i] - mean,0));
            }

            // Get maxium lag
            int maxlag = y.Count() - 1;
            // Perform the auto- and cross-correlations.

            // Doubling length of y
            int y_length = 2 * y.Count();
            for (int i = y.Count(); i < y_length; i++)
            {
                y.Add(new Complex(0, 0));
            }



            Complex[] y_complex = y.ToArray();

            MathNet.Numerics.IntegralTransforms.Fourier.Forward(y_complex); // In the matlab script, a transform length of twice the array length was used. I dont know how to do this and dont know if it is necessary
  

            // correlation list
            List<Complex> Cr = new List<Complex>();
            for (int i = 0; i < y_complex.Length; i++)
            {
                Cr.Add(new Complex((y_complex[i].Magnitude * y_complex[i].Magnitude),0)); // Abs(X).^2
            }
            List<double> Cr_double = new List<double>();
            for (int i = 0; i < Cr.Count(); i++)
            {
                Cr_double.Add(Cr[i].Real);
            }
            save_to_csv(Cr_double, @"C:\github\PCG_Segmentation\Cr_double.csv");

            Complex[] Cr_array = Cr.ToArray();

            MathNet.Numerics.IntegralTransforms.Fourier.Inverse(Cr_array); // In the matlab script, a transform length of twice the array length was used. I dont know how to do this and dont know if it is necessary

            List<double> Cr_ifft_double = new List<double>();
            for (int i = 0; i < Cr.Count(); i++)
            {
                Cr_ifft_double.Add(Cr_array[i].Real);
            }
            save_to_csv(Cr_ifft_double, @"C:\github\PCG_Segmentation\Cr_ifft_double.csv");
            // Dividing by max to normalise it (IS THIS WHAT IM MEANT TO DO???)
            double max_Cr_ifft = Cr_ifft_double.Max();
            List<double> Cr_ifft_shifted = new List<double>(); // Shifing Cr
            for (int i =Cr_ifft_double.Count()/2; i<Cr_ifft_double.Count(); i++)
            {
                Cr_ifft_shifted.Add(Cr_ifft_double[i]/max_Cr_ifft);
            }
            for (int i = 0; i < Cr_ifft_double.Count() / 2; i++)
            {
                Cr_ifft_shifted.Add(Cr_ifft_double[i]/ max_Cr_ifft);
            }
            save_to_csv(Cr_ifft_shifted, @"C:\github\PCG_Segmentation\Cr_ifft_shifted.csv");

            // Get signal autocorrelation 
            List<double> auto_correlation = new List<double>();
            for (int i = Cr_ifft_shifted.Count()/2; i < Cr_ifft_shifted.Count(); i++)
            {
                auto_correlation.Add(Cr_ifft_shifted[i]);
            }

            // Low pass filter the autocorrelation

            if (Fs >= 15)
            {
                LowpassFilterButterworthImplementation filter_low = new LowpassFilterButterworthImplementation(15, 2, Fs);
                for (int i = 0; i < auto_correlation.Count(); i++)
                {
                    signal_autocorrelation_filtered.Add(filter_low.compute(auto_correlation[i]));
                }
                save_to_csv(signal_autocorrelation_filtered, @"C:\github\PCG_Segmentation\signal_autocorrelation_filtered.csv");

            }

        }
        */
        /*
        public static double get_hr_peak_autocorrelation(List<double> signal_autocorrelation, double max_HR, double min_HR, int Fs)
        {
            int min_index = Convert.ToInt32(Math.Round((60 / max_HR) * Fs));
            int max_index = Convert.ToInt32(Math.Round((60 / min_HR) * Fs));

            double max_val = 0;
            int index = 0;
            for (int i = min_index; i < max_index; i++)
            {
                if (signal_autocorrelation[i] > max_val)
                {
                    max_val = signal_autocorrelation[i];
                    index = i;
                }
            }

            double HR = 60.0 / (Convert.ToDouble(index) / Convert.ToDouble(Fs));
            return HR;
        }

        public static void get_systolicTimeInterval(double systolicTimeInterval, List<double> signal_autocorrelation, double heartRate, double max_HR, int Fs)
        {
            int max_sys_duration = Convert.ToInt32((60.0 / heartRate) * Convert.ToDouble(Fs) / 2.0 );
            double sysHR = max_HR + 30.0;

            int min_sys_duration = Convert.ToInt32((60.0 / sysHR) * Convert.ToDouble(Fs) / 2.0);

            double max_val = 0;
            double index = 0;
            for (int i = min_sys_duration; i < max_sys_duration; i++)
            {
                double temp = signal_autocorrelation[i]; 
                if (signal_autocorrelation[i] > max_val)
                {
                    max_val = signal_autocorrelation[i];
                    index = Convert.ToDouble(i);
                }
            }

            systolicTimeInterval = (index) / Convert.ToDouble(Fs);
        }

        */
        /*
        public static double[][] Multiply(double[][] matrix1, double[][] matrix2)
        {
            // cahing matrix lengths for better performance  
            var matrix1Rows = matrix1.Length;
            var matrix1Cols = matrix1[0].Length;
            var matrix2Rows = matrix2.Length;
            var matrix2Cols = matrix2[0].Length;

            // checking if product is defined  
            if (matrix1Cols != matrix2Rows)
                throw new InvalidOperationException
                  ("Product is undefined. n columns of first matrix must equal to n rows of second matrix");

            // creating the final product matrix  
            //double[][] product = new double[matrix1Rows, matrix2Cols];
            double[][] product = new double[matrix1Rows][];
            for (int i = 0; i < matrix1Rows; i++)
            {
                product[i] = new double[matrix2Cols];
            }
            // looping through matrix 1 rows  
            for (int matrix1_row = 0; matrix1_row < matrix1Rows; matrix1_row++)
            {
                // for each matrix 1 row, loop through matrix 2 columns  
                for (int matrix2_col = 0; matrix2_col < matrix2Cols; matrix2_col++)
                {
                    // loop through matrix 1 columns to calculate the dot product  
                    for (int matrix1_col = 0; matrix1_col < matrix1Cols; matrix1_col++)
                    {
                        product[matrix1_row][matrix2_col] +=
                          matrix1[matrix1_row][matrix1_col] *
                          matrix2[matrix1_col][matrix2_col];
                    }
                }
            }

            return product;
        }
        */


        // This is a function that computes the Predict values for a nominal or ordinal multinomial regression model.
        // based on the matlab function mnrval

        /*
        public static double[][] mnrval(double[] beta, double[][] x)
        {
            // Validate the size of beta and compute the linear predictors

            // matlab: eta = repmat(beta(1:(k-1))',n,1) + repmat(x*beta(k:pstar),1,k-1);

            int pstar = beta.Length;
            int p = x.Length;
            int k = pstar - p;

            // find repmat(beta(1:(k-1))',n,1)

            double[][] repmat1 = new double[x[0].Length][];
            for (int i = 0; i < x[0].Length; i++)
            {
                repmat1[i] = new double[k];
            }

            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < x[0].Length; j++)
                {
                    repmat1[j][i] = beta[i];
                }
            }

            double[][] beta_k = new double[beta.Length - k][];
            for (int i = 0; i < beta.Length-k; i++)
            {
                beta_k[i] = new double[k];
            }
            for (int i = 0; i < k; i++)
            {
                for (int j = k; j < beta.Length; j++)
                {
                    beta_k[j-1][i] = beta[j];
                }
            }

            double[][] x_transform = Transform(x);
            double[][] x_times_beta = Multiply(x_transform, beta_k);
            double[][] eta = new double[x_times_beta.Length][];
            for (int i =0; i < eta.Length; i++)
            {
                eta[i] = new double[x_times_beta[i].Length];
            }
            for (int i = 0; i < eta.Length; i++)
            {
                for (int j = 0; j < eta[i].Length; j++)
                {
                    eta[i][j] = x_times_beta[i][j] + repmat1[i][j];
                }
            }


            // Using only case nominal

            // finding max of each row
            double[] maxeta = new double[eta.Length];

            for (int i = 0; i < eta.Length; i++)
            {
                maxeta[i] = double.MinValue;
                for (int j = 0; j < eta[i].Length; j++)
                {
                    if (eta[i][j] > maxeta[i])
                    {
                        maxeta[i] = eta[i][j];
                    }
                }
            }

            // rescale so maximum probability is 1
            double[][] pi = new double[maxeta.Length][];
            for (int i = 0; i < maxeta.Length; i++)
            {
                // rescale so max probability is 1
                pi[i] = new double[2];
                pi[i][0] = Math.Exp(eta[i][0] - maxeta[i]);
                pi[i][1] = Math.Exp(- maxeta[i]);

                // renormalise for real probabilities
                double temp = pi[i][0];
                pi[i][0] = pi[i][0] / (temp + pi[i][1]);
                pi[i][1] = pi[i][1] / (temp + pi[i][1]);
            }


            return pi;

        }

        public static double mvnpdf(double X, double Mu, double  Sigma)
        {
            double X0 = X;
            double xRinv = X;
            double logSqrtDetSigma = 0;
            double d = 1;

            // Get vector mean
            X0 = X - Mu;


            // Hardcoding R

            double R = Sigma;
            xRinv = X0 / R;
            logSqrtDetSigma = Math.Log(R);

            double quadform = Math.Pow(xRinv, 2);
            double y = Math.Exp(-0.5 * quadform - logSqrtDetSigma - d * Math.Log(2.0 * Math.PI) / 2.0);
            return y;



        }

        public static double mvnpdf(double[] X, double[] Mu, double[][] Sigma)
        {

            double logSqrtDetSigma = 0;
            double[][] X0 = new double[1][];
            X0[0] = new double[X.Length];

            // Get vector mean
            for (int i = 0; i < X0[0].Length; i++)
            {
                X0[0][i] = X[i] - Mu[i];
            }

            // Make sure Sigma is a valid covariance matrix

            // Matlab calculates the cholcov: [R,err] = cholcov_copied(Sigma,0);
            // We will hardcode this as it uses Sigma which is always the same




            double[][] R = new double[4][];
            R[0] = new double[] { 0.9992, 0.9043, 0.7432, 0.7075 };
            R[1] = new double[] { 0, 0.4251, 0.4196, 0.3632 };
            R[2] = new double[] { 0, 0, 0.5197, 0.0551 };
            R[3] = new double[] { 0, 0, 0, 0.6025 };



            // Create array of standardized data, and compute log(sqrt(det(Sigma)))

            double[][] xRinv = Matrix.Divide(X0, R);

            // logSqrtDetSigma = sum(log(diag(R)));
            for (int i = 0; i < R.Length; i++)
            {
                logSqrtDetSigma += Math.Log(R[i][i]);
            }

            // % The quadratic form is the inner products of the standardized data
            // quadform = sum(xRinv.^ 2, 2);
            double quadform = 0;
            for (int i = 0; i < xRinv.Length; i++)
            {
                for (int j = 0; j < xRinv[i].Length; j++)
                {
                    quadform += Math.Pow(xRinv[i][j], 2);
                }
            }

            // y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
            int d = X.Length;
            double y = Math.Exp(-0.5 * quadform - logSqrtDetSigma - Convert.ToDouble(d) * Math.Log(2.0 * Math.PI) / 2.0);
            return y;
        }

        public static double[] getColumn(double[][] input, int desired_column)
        {
            double[] result = new double[input.Length];
            for (int i = 0; i < input.Length; i++)
            {
                result[i] = input[i][desired_column];
            }

            return result;
        }

        public static double[][] Transform(double[][] input)
        {
            double[][] result = new double[input[0].Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[input.Length];
            }

            int num_rows = result.Length;
            int num_columns = input.Length;

            for (int i = 0; i < (num_rows * num_columns); i++)
            {
                result[i / num_columns][i % num_columns] = input[i % num_columns][i / num_columns];
            }

            return result;

        }
        */

        public static List<double> smooth_data(List<double> input, int num_smooth)
        {
            List<double> output = new List<double>();
            

            for (int i = 0; i < num_smooth / 2; i++)
            {
                output.Add(0);
            }
            for (int i = num_smooth/2; i < input.Count() - num_smooth/2; i++)
            {
                double average = 0;
                for (int j = i - num_smooth /2; j <= i + num_smooth / 2; j++)
                {
                    average += input[j];
                }
                average /= Convert.ToDouble(num_smooth);
                output.Add(average);
            }

            for (int i = input.Count() - num_smooth/2; i < input.Count(); i++)
            {
                output.Add(0);
            }

            return output;

        }


        public static (List<double>, List<double>s) Noisy_Peak_Finding(List<double> input, int num_smooth, double threshold)
        {
            
            //input = smooth_data(input, 10);

            //PCG_methods.save_to_csv(smoothed_signal, @"C:\github\PCG_Segmentation\smoothed_signal.csv");
            //PCG_methods.save_to_csv(input, @"C:\github\PCG_Segmentation\input_to_.csv");
            //Viterbi_Springer.PlotData(smoothed_signal, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\smoothed_signal.png");


            List<double> peakIndices = new List<double>();
            List<double> peakValues = new List<double>();
            

            

            double peak_Value = double.MinValue;
            int peak_Index = int.MinValue;

            for (int i = 0; i < input.Count(); i++)
            {
                double value = input[i];
                if (value > threshold)
                {
                    if ((peak_Value == double.MinValue) || (value > peak_Value))
                    {
                        peak_Index = i;
                        peak_Value = value;
                    }
                }
                else if ((value < threshold) & peak_Index != int.MinValue)
                {
                    peakIndices.Add(peak_Index);
                    peakValues.Add(peak_Value);
                    peak_Value = double.MinValue;
                    peak_Index = int.MinValue;
                }

            
            }

            //if (peak_Index != int.MinValue)
            //{
            //    peakIndices.Add(peak_Index);
            //}
            

            return (peakIndices, peakValues);

        }

        public static (List<double>, List<double>) peak_post_processing(List<double> input, List<double> input_Y, int Fs)
        {
            List<double> S1_timings = new List<double>();
            List<double> S1_Amplitudes = new List<double>();
            List<double> S2_timings = new List<double>();
            List<double> S2_Amplitudes = new List<double>();

            int S_flag = 0;
            int uncertain_flag = 1;
            int certain_amount = 3;

            double min_difference = 0.1;



            // Check that all samples are valid

            for (int i = 0; i < input.Count() - 1; i++)
            {
                double min_diff = (min_difference) * Convert.ToDouble(Fs);
                double diff = Convert.ToDouble(input[i + 1] - input[i]);
                if ((diff) < min_diff)
                {
                    input.RemoveAt(i + 1);
                    input_Y.RemoveAt(i + 1);
                }
            }
            /*
            for (int i = 0; i < input.Count() - 1; i++)
            {
                if (input_Y[i] > input_Y[i + 1])
                {
                    S2_timings.Add(input[i]);
                    S2_Amplitudes.Add(input_Y[i]);
                    S_flag = 1;
                }
                else
                {
                    S1_timings.Add(input[i]);
                    S1_Amplitudes.Add(input_Y[i]);
                    S_flag = 2;
                }
            }
            if (S_flag == 1)
            {
                S1_timings.Add(input[input.Count() - 1]);
                S1_Amplitudes.Add(input_Y[input.Count() - 1]);
            }
            else
            {
                S2_timings.Add(input[input.Count() - 1]);
                
                S2_Amplitudes.Add(input_Y[input.Count() - 1]);
            }
            */

            return (input, input_Y);
        }
        public static List<double> samples2Timings(List<double> input, int Fs)
        {
            List<double> timings = new List<double>();
            for (int i = 0; i < input.Count(); i++)
            {
                timings.Add((input[i]) / Fs);
            }
            return timings;
        }


        public static List<double> crop_data(List<double> input, int input_Fs, int cut_num, int cut_Fs)
        {
            List<double> output = new List<double>();
            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i]);
            }
            output.RemoveRange(0, (cut_num * input_Fs / cut_Fs) / 2);
            output.RemoveRange(output.Count() - (cut_num * input_Fs / cut_Fs) / 2, (cut_num * input_Fs / cut_Fs) / 2);


            return output;
        }
    }






}
