using System;
using System.Collections.Generic;
using System.Linq;
using Spectrogram;
using ScottPlot;
using Accord.Math;


//using WaveletStudio;
//using WaveletStudio.Blocks;

namespace PCG_segmentation
{
    internal class Viterbi_Springer
    {

        public static void PlotData(List<double> data, double samplingFrequency, string ylabel, string filepath)
        {
            var plt = new Plot(800, 600);
            plt.AddSignal(data.ToArray(), samplingFrequency);
            plt.YLabel(ylabel);
            plt.Margins(0);
            plt.SaveFig(filepath);
        }

        public static void MultiPlot(List<double> data1, List<double> data2X, List<double> data2Y, string filepath)
        {
            var plt = new Plot(800, 600);
            var mySignalPlot1 = plt.AddSignal(data1.ToArray());
            mySignalPlot1.YAxisIndex = 0;
            mySignalPlot1.XAxisIndex = 0;

            var mySignalPlot2 = plt.PlotScatter(data2X.ToArray(), data2Y.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            plt.SaveFig(filepath);
        }

        public static void MultiPlot(List<double> data1, List<double> data2, string filepath)
        {
            var plt = new Plot(800, 600);
            var mySignalPlot1 = plt.AddSignal(data1.ToArray());
            mySignalPlot1.YAxisIndex = 0;
            mySignalPlot1.XAxisIndex = 0;

            var mySignalPlot2 = plt.AddSignal(data2.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            plt.SaveFig(filepath);
        }

        public static void TriplePlot(List<double> data1, List<double> data2X, List<double> data2Y, List<double> data3X, List<double>data3Y, string filepath)
        {
            var plt = new Plot(800, 600);
            var mySignalPlot1 = plt.AddSignal(data1.ToArray());
            mySignalPlot1.YAxisIndex = 0;
            mySignalPlot1.XAxisIndex = 0;

            var mySignalPlot2 = plt.PlotScatter(data2X.ToArray(), data2Y.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            var mySignalPlot3 = plt.PlotScatter(data3X.ToArray(), data3Y.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            plt.SaveFig(filepath);
        }
        static void resample(List<double> output, List<double> input, int current_Fs, int new_FS)
        {
            for (int i = 0; i < current_Fs; i++)
            {
                if (i % (current_Fs / new_FS) == 0)
                {
                    output.Add(input[i]);

                }
            }
        }

        public static List<double> Resample(List<double> inputList, int outputLength)
        {
            List<double> outputList = new List<double>(outputLength);

            for (int i = 0; i < outputLength; i++)
            {
                double index = i * (inputList.Count - 1) / (double)(outputLength - 1);

                int lowerIndex = (int)Math.Floor(index);
                int upperIndex = (int)Math.Ceiling(index);

                double lowerValue = inputList[lowerIndex];
                double upperValue = inputList[upperIndex];

                double interpolatedValue = lowerValue + (upperValue - lowerValue) * (index - lowerIndex);

                outputList.Add(interpolatedValue);
            }

            return outputList;
        }


        /*
        public static (double[][], double, double, double, double, double, double, double, double) get_duration_distributions_modified(double heartrate, double systolic_time)
        {
            // Loading default springer options

            double audio_Fs = 1000.0;
            double audio_segmentation_Fs = 50;
            double segmentation_tolerance = 0.1;
            int use_mex = 1;
            int include_wavelt_feature = 1;


            double mean_factor = 1;
            double std_factor = 1;
            int mean_s1 = Convert.ToInt32(0.122 / mean_factor * audio_segmentation_Fs);
            int std_s1 = Convert.ToInt32(0.022 * std_factor * audio_segmentation_Fs);

            int mean_s2 = Convert.ToInt32(0.094 / mean_factor * audio_segmentation_Fs);
            int std_s2 = Convert.ToInt32(0.022 * std_factor * audio_segmentation_Fs);

            int mean_systole = Convert.ToInt32(systolic_time * audio_segmentation_Fs) - mean_s1;
            double std_systole = (25.0 / 1000.0) * audio_segmentation_Fs;

            double mean_diastole = ((60 / heartrate) - systolic_time - 0.094) * audio_segmentation_Fs;
            double std_diastole = 0.07 * mean_diastole + (6 / 1000) * audio_segmentation_Fs;

            // %% Cell array for the mean and covariance of the duration distributions:

            int num_rows = 4;
            int num_columns = 2;
            double[][] d_distribution = new double[num_rows][];
            for (int i = 0; i < num_rows; i++)
            {
                d_distribution[i] = new double[num_columns];
            }

            d_distribution[0][0] = mean_s1;
            d_distribution[0][1] = Math.Pow(std_s1, 2);
            d_distribution[1][0] = mean_systole;
            d_distribution[1][1] = Math.Pow(std_systole, 2);
            d_distribution[2][0] = mean_s2;
            d_distribution[2][1] = Math.Pow(std_s2, 2);
            d_distribution[3][0] = mean_diastole;
            d_distribution[3][1] = Math.Pow(std_diastole, 2);

            // %Min systole and diastole times
            double min_systole = mean_systole - 3 * (std_systole + std_s1);
            double max_systole = mean_systole + 3 * (mean_systole + std_s1);

            double min_diastole = mean_diastole - 3 * std_diastole;
            double max_diastole = mean_diastole + 3 * std_diastole;

            // %Setting the Min and Max values for the S1 and S2 sounds:
            // % If the minimum lengths are less than a 50th of the sampling frequency, set
            // % to a 50th of the sampling frequency:

            double min_s1 = (mean_s1 - 3 * (std_s1));

            if (min_s1 < audio_segmentation_Fs / 50.0)
            {
                min_s1 = audio_segmentation_Fs / 50.0;
            }

            double min_s2 = (mean_s2 - 3 * (std_s2));
            if (min_s2 < (audio_segmentation_Fs / 50))
            {
                min_s2 = (audio_segmentation_Fs / 50);
            }

            double max_s1 = (mean_s1 + 3 * (std_s1));
            double max_s2 = (mean_s2 + 3 * (std_s2));

            return (d_distribution, max_s1, min_s1, max_s2, min_s2, max_systole, min_systole, max_diastole, min_diastole);


        }

        */

        public static (List<double>, int) get_PSD_feature_Springer_HMM(List<double> data, int sampling_frequency, int frequency_limit_lower, int frequency_limit_upper)
        {



            //PlotData(data, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PSD_input_data.png");
            // Formuila for number of ffts size: fix((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP))

            var sg = new SpectrogramGenerator(sampleRate: sampling_frequency, fftSize: 128, stepSize: sampling_frequency / 80); // I am not sure what an appropropriate step size is
            //var sg = new SpectrogramGenerator(sampleRate: sampling_frequency, fftSize: 2048, stepSize: 1);
            double[] window = { 0.0800, 0.0957, 0.1416, 0.2147, 0.3100, 0.4209, 0.5400, 0.6591, 0.7700, 0.8653, 0.9384, 0.9843, 1.0000, 0.9843, 0.9384, 0.8653, 0.7700, 0.6591, 0.5400, 0.4209, 0.3100, 0.2147, 0.1416, 0.0957, 0.0800 }; // change this to generating the hanning window 
            
            sg.SetWindow(window);
            double[] data_array = data.ToArray();
            //PCG_methods.save_to_csv(data, @"C:\github\PCG_Segmentation\THIS_SHOULD_BE_NORMAL.csv");
            sg.Add(data_array);

            //sg.SaveImage(@"C:\github\PCG_Segmentation\PCG_segmentation\hal.png");

            List<double[]> FFTs = sg.GetFFTs();

            for (int i = 0; i < FFTs.Count(); i++)
            {
                for (int j = 0; j < FFTs[i].Length; j++)
                {
                    FFTs[i][j] = Math.Pow(Math.Abs(FFTs[i][j]), 2) / (9.544) * 1000;
                }
            }

            // Find index with frequency = 40. HzPerPx = 0.97 -> 0.97 * x = 40 -> x = 40/0.97

            int lower_limit_position = Convert.ToInt32(40 / sg.HzPerPx);
            int upper_limit_position = Convert.ToInt32(60 / sg.HzPerPx);

            List<double> psd = new List<double>();
            for (int i = 0; i < FFTs.Count(); i++)
            {
                psd.Add(0.0);
                for (int j = lower_limit_position; j < upper_limit_position; j++)
                {
                    psd[i] = psd[i] + FFTs[i][j];
                }
                psd[i] = psd[i] / (upper_limit_position - lower_limit_position);
            }


            /*
            // Take absolute value of each point of FFT. Square it. divide by sampling frequency and squared sum of the window i.e. divide by (9.544 * 1000)
            // sum across frequencies 40 to 60 hz. that is the power spectral envelope.
            var plt = new ScottPlot.Plot(800, 600);
            plt.AddSignal(psd.ToArray(), sampling_frequency);
            plt.YLabel("Amplitude");
            plt.Margins(0);
            plt.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\psd.png");

            data_array = data.ToArray();
            plt = new ScottPlot.Plot(800, 600);
            plt.AddSignal(data_array, sampling_frequency / 1000.0);
            plt.YLabel("Amplitude");
            plt.Margins(0);
            plt.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\psd_input_data.png");
            */
            // Data is cut when performing the PSD. here I am concatenating 0s to each end so it remains the same length.

            // Amount of data cut when Fs = 1000
            int cut_data = data.Count() - (psd.Count() * sg.StepSize);
            // Amount of data cut when Fs = 1000/12
            int toAppend = Convert.ToInt32(Convert.ToDouble(data.Count()) / 12.0) - psd.Count();
            for (int i = 0; i < toAppend / 2; i++)
            {
                psd.Insert(0, 0.0);
                psd.Add(0.0);
            }
            

            return (psd, cut_data);

        }
        /*
        static public void viterbiDecodePCG_Springer_modified(List<List<double>> observation_sequence, double[] pi_vector, double[][] b_matrix, double[] total_obs_distribution1, double[][] total_obs_distribution2, double heartrate, double systolic_time, int Fs)
        {
            // Using default springer_options

            int audio_Fs = 1000;
            int audio_segmentation_Fs = 50;
            double segmentation_tolerance = 0.1;
            int use_mex = 1;

            // we do not have a wavelet feature. so we will set this to 0

            int include_wavelet_feature = 0;

            int T = observation_sequence[0].Count();
            int N = 4; // The number of states

            // Setting the maximum duration of a single state. This is set to an entire
            // heart cycle:

            int max_duration_D = Convert.ToInt32(1.0 * (60.0 / heartrate) * Convert.ToDouble(audio_segmentation_Fs));

            // Initialising the variables that are needed to find the optimal state path along
            // the observation sequence.
            // delta_t(j), as defined on page 264 of Rabiner, is the best score(highest
            // probability) along a single path, at time t, which accounts for the first
            // t observations and ends in State s_j.In this case, the length of the
            // matrix is extended by max_duration_D samples, in order to allow the use
            // of the extended Viterbi algortithm:

            
            double[][] delta = new double[max_duration_D][];
            for (int i =0; i < max_duration_D; i++)
            {
                delta[i] = new double[N];
                for(int j = 0; j< N; j++)
                {
                    delta[i][j] = double.MinValue;
                }
            }


            // The argument that maximises the transition between states(this is
            // basically the previous state that had the highest transition probability
            // to the current state) is tracked using the psi variable.

            List<List<double>> psi_duration = new List<List<double>>();

            for (int i = 0; i < N; i++)
            {
                psi_duration.Add(new List<double>());
                for (int j = 0; j < T + max_duration_D; j++)
                {
                    psi_duration[i].Add(0.0);
                }
            }

            // Setting up observation probs

            List<List<double>> observation_probs = new List<List<double>>();

            for (int i = 0; i < N; i++)
            {
                observation_probs.Add(new List<double>());
                for (int j = 0; j < T; j++)
                {
                    observation_probs[i].Add(0.0);
                }
            }


            // double[][] observation_probs_array = observation_probs.Select(list => list.ToArray()).ToArray();
            double[][] observation_probs_array = new double[observation_probs.Count()][];

            for (int i =0; i < N; i++)
            {
                observation_probs_array[i] = observation_probs[i].ToArray();
            }

            // Converting observation_sequence to array

            double[][] observation_sequence_array = new double[observation_sequence.Count()][];
            for (int i = 0; i < observation_sequence.Count(); i++)
            {
                observation_sequence_array[i] = observation_sequence[i].ToArray();
            }
            
            for (int i = 0; i < N; i++)
            {
                double[][] pihat = PCG_methods.mnrval(b_matrix[i], observation_sequence_array);

                for (int j = 0; j < T; j++)
                {
                    double[] curr_observation_sequence = new double[observation_sequence.Count()];
                    for (int k = 0; k < curr_observation_sequence.Length; k++)
                    {
                        curr_observation_sequence[i] = observation_sequence[k][j];
                    }
                    
                    double Po_correction = PCG_methods.mvnpdf(curr_observation_sequence, total_obs_distribution1, total_obs_distribution2);

                    //% When saving the coefficients from the logistic
                    //% regression, it orders them P(class 1) then P(class 2). When
                    //%training, I label the classes as 0 and 1, so the
                    //%correct probability would be pihat(2).

                    observation_probs[i][j] = pihat[j][1] * Po_correction / pi_vector[i];
                }


            }

            // %% Setting up state duration probabilities, using Gaussian distributions:
            (double[][] d_distributions, double max_S1, double min_S1, double max_S2, double min_S2, double max_systole, double min_systole, double max_diastole, double min_diastole) = get_duration_distributions_modified(heartrate, systolic_time);

            double[][] duration_probs = new double[N][];
            for (int i = 0; i < N; i++)
            {
                duration_probs[i] = new double[3 * audio_segmentation_Fs];
                for(int j =0; j< audio_segmentation_Fs; j++)
                {
                    duration_probs[i][j] = 0;
                }
                
            }

            double[] duration_sum = new double[N];
            for (int i = 0; i < N; i ++)
            {
                duration_sum[i] = 0;
            }
            for (int state_j = 0; state_j < N; state_j++)
            {
                for (int d = 0; d < max_duration_D; d++)
                {
                    if (state_j == 0)
                    {
                        duration_probs[state_j][d] = PCG_methods.mvnpdf(Convert.ToDouble(d), d_distributions[state_j][0], d_distributions[state_j][1]);
                        if (d < min_S1 || d > max_S1)
                        {
                            duration_probs[state_j][d] = double.NegativeInfinity;
                        }
                    }
                    else if (state_j == 1)
                    {
                        duration_probs[state_j][d] = PCG_methods.mvnpdf(Convert.ToDouble(d), d_distributions[state_j][0], d_distributions[state_j][1]);
                        if (d < min_systole || d > max_systole)
                        {
                            duration_probs[state_j][d] = double.NegativeInfinity;
                        }
                    }
                    else if (state_j == 2)
                    {
                        duration_probs[state_j][d] = PCG_methods.mvnpdf(Convert.ToDouble(d), d_distributions[state_j][0], d_distributions[state_j][1]);
                        if (d < min_S2 || d > max_S2)
                        {
                            duration_probs[state_j][d] = double.NegativeInfinity;
                        }
                    }
                    else if (state_j == 3)
                    {
                        duration_probs[state_j][d] = PCG_methods.mvnpdf(Convert.ToDouble(d), d_distributions[state_j][0], d_distributions[state_j][1]);
                        if (d < min_diastole || d > max_diastole)
                        {
                            duration_probs[state_j][d] = double.NegativeInfinity;
                        }
                    }

                    duration_sum[state_j] += duration_probs[state_j][d];
                }
                
            }

            // there is this instruction: 

            //if (length(duration_probs) > 3 * Fs)
            //    duration_probs(:, (3 * Fs + 1):end) = [];
            //end

            // but it seems to serve no purpose

            // PERFORMING THE ACTUAL VITERBI RECURSION

            double[] qt = new double[delta[0].Length];
            double[][] test = new double[0][];
            

            // Initialisation step
            

        }
        */
     }
}