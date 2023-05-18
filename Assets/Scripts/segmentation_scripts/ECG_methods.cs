using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using DSP;
using FftSharp;
using ScottPlot;

namespace PCG_segmentation
{
    internal class ECG_methods
    {
        public static List<double> ecg_clean(List<double> ecg_signal, int sampling_rate)
        {
            List<double> cleaned_data = new List<double>();

            //LowpassFilterButterworthImplementation filter_low = new LowpassFilterButterworthImplementation(400, 2, Fs_new);

            // Remove baseline wander
            HighpassFilterButterworthImplementation filter_high = new HighpassFilterButterworthImplementation(0.5, 2, sampling_rate);
            for (int i = 0; i < ecg_signal.Count(); i++)
            {
                cleaned_data.Add(filter_high.compute(ecg_signal[i]));
            }

            // Remove high frequency components

            LowpassFilterButterworthImplementation filter_Low = new LowpassFilterButterworthImplementation(150, 3, sampling_rate);
            //for (int i = 0; i < ecg_signal.Count(); i++)
            //{
            //    cleaned_data.Add(filter_Low.compute(ecg_signal[i]));
            //}

            for (int i = 0; i < ecg_signal.Count(); i++)
            {
                cleaned_data[i] = (filter_Low.compute(cleaned_data[i]));
            }

            return cleaned_data;
        }

        public static List<double> median_filter(List<double> input, int width)
        {
            List<double> output = new List<double>();


            for (int i = 0; i < width / 2; i++)
            {
                output.Add(0);
            }
            for (int i = width / 2; i < input.Count() - width / 2; i++)
            {
                List<double> median_list = new List<double>();

                for (int j = i - width / 2; j <= i + width / 2; j++)
                {
                    median_list.Add(input[j]);
                }
                median_list.Sort();
                output.Add(median_list[width / 2]);
            }

            for (int i = input.Count() - width / 2; i < input.Count(); i++)
            {
                output.Add(0);
            }

            return output;
        }

        public static (List<double>, List<double>) ECG_peak_detection(List<double> input, double threshold)
        {
            List<double> peakIndices = new List<double>();
            List<double> peakValues = new List<double>();

            double peak_Value = double.MinValue;
            int peak_Index = int.MinValue;
            // Finding based on threshold
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

            return (peakIndices, peakValues);
        }

        public static List<double> Remove_R_Peaks(List<double> input, List<double> R_indices, int window)
        {
            bool positive_flag = false;
            int k;

            List<double> output = new List<double>();
            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i]);
            }
            for (int i = 0; i < R_indices.Count(); i++)
            {
                if (Convert.ToInt32(R_indices[i]) - window  < 0)
                {
                    for (int j = 0; j < Convert.ToInt32(R_indices[i]) + window; j++)
                    {
                        output[j] = 0.0;
                    }
                }
                else if (Convert.ToInt32(R_indices[i]) + window > input.Count())
                {
                    for (int j = Convert.ToInt32(R_indices[i]) - window ; j < input.Count(); j++)
                    {
                        output[j] = 0.0;
                    }

                    
                }
                else
                {
                    for (int j = Convert.ToInt32(R_indices[i]) - window ; j < Convert.ToInt32(R_indices[i]) + window; j++)
                    {
                        output[j] = 0.0;
                    }
                    
                    if (output[Convert.ToInt32(R_indices[i]) - window - 1] > 0)
                    {
                        k = Convert.ToInt32(R_indices[i]) - window - 1;
                        positive_flag = true;
                        while (positive_flag & k > 0)
                        {
                            output[k] = 0.0;
                            k--;
                            if (output[k] < 0 | output[k] > 0.5)
                            {
                                positive_flag = false;
                            }
                        }

                    }
                    
                }



            }



            return output;

        }

        static public (List<double>, List<double>, List<double>, List<double>) post_processing(List<double> R_timings, List<double> T_timings, List<double> S_timings, bool figures)
        {
            // The order we want is P_wave, R_peaks, S1, T_peaks, S2
            // For now I am not using P_wave
            List<double> S1_timings = new List<double>();
            List<double> S2_timings = new List<double>();
            List<double> T_timings_out = new List<double>();
            double min_HR = 40.0; // BPM
            double max_diff = min_HR / 60.0;
            double max_HR = 200.0; // BPM
            int num_missed_S = 0;
            int num_missed_T = 0;
            double tolerance = 0.2;



            double expected_S1_R_diff = 0.075;
            double expected_S2_R_diff = 0.33;
            double expected_T_R_diff = 0.199;
            // Remove initial timings so that R_peak is the first timing

            if (T_timings[0] < R_timings[0])
            {
                T_timings.RemoveAt(0);
            }
            if (S_timings[0] < R_timings[0])
            {
                S_timings.RemoveAt(0);
            }

            // Iterate through List and ensure order is correct

            // Ensure that no samples have been missed. 

            // there should be a MINIMUM JUMP and a MAXIMUM JUMP
            R_timings.Add(2 * R_timings[R_timings.Count() - 1] - R_timings[R_timings.Count() - 2]);
            for (int i = 0; i < R_timings.Count()-1 ; i++)
            {

                int num_T_peaks = 0;

                int T_j = 0;


                if (T_timings.Count() == 0)
                {
                    num_T_peaks = 0;
                    T_j = 0;
                }
                else
                {

                    while (T_timings[T_j] < R_timings[i + 1])
                    {
                        num_T_peaks++;
                        T_j++;
                        if (T_j >= T_timings.Count())
                        {
                            break;
                        }
                    }
                }


                if (num_T_peaks == 0)
                {
                    T_timings_out.Add(R_timings[i] + expected_T_R_diff);
                    num_missed_T++;
                }
                else if (num_T_peaks == 1)
                {
                    if (T_timings[0] < (R_timings[i] + expected_T_R_diff + tolerance) & T_timings[0] > R_timings[i] + expected_T_R_diff - tolerance)
                    {
                        T_timings_out.Add(T_timings[0]);
                    }
                    else
                    {
                        T_timings_out.Add(R_timings[i] + expected_T_R_diff);
                        num_missed_T++;

                    }
                    T_timings.RemoveAt(0);
                }
                else
                {
                    T_timings_out.Add(R_timings[i] + expected_T_R_diff);
                    num_missed_T++;
                    for (int k = 0; k < num_T_peaks; k++)
                    {
                        T_timings.RemoveAt(0);
                    }
                }

                int num_peaks = 0;
                int j = 0;

                if (S_timings.Count() == 0)
                {
                    num_peaks = 0;
                    j = 0;
                }
                else
                {
                    while (S_timings[j] < R_timings[i + 1])
                    {
                        num_peaks++;
                        j++;
                        if (j >= S_timings.Count())
                        {
                            break;
                        }
                    }
                }

                if (num_peaks == 0)
                {
                    num_missed_S++;
                    S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                    S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                }
                else if (num_peaks == 1)
                {
                    if (S_timings[0] < (R_timings[i] + expected_S1_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S1_R_diff - tolerance)
                    {
                        num_missed_S++;
                        S1_timings.Add(S_timings[0]);
                        S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    }
                    else if (S_timings[0] < (R_timings[i] + expected_S2_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S2_R_diff - tolerance)
                    {
                        num_missed_S++;
                        S2_timings.Add(S_timings[0]);
                        S1_timings.Add(R_timings[i] + expected_S1_R_diff);

                    }
                    S_timings.RemoveAt(0);
                }
                else if (num_peaks == 2)
                {
                    if (S_timings[0] < (R_timings[i] + expected_S1_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S1_R_diff - tolerance)
                    {
                        S1_timings.Add(S_timings[0]);
                        if (S_timings[1] < (R_timings[i] + expected_S2_R_diff + tolerance) & S_timings[1] > R_timings[i] + expected_S2_R_diff - tolerance)
                        {
                            S2_timings.Add(S_timings[1]);
                        }
                        else
                        {
                            num_missed_S++;
                            S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                        }
                    }
                    else if (S_timings[0] < (R_timings[i] + expected_S2_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S2_R_diff - tolerance)
                    {
                        S2_timings.Add(S_timings[0]);
                        if (S_timings[1] < (R_timings[i] + expected_S1_R_diff + tolerance) & S_timings[1] > R_timings[i] + expected_S1_R_diff - tolerance)
                        {
                            S1_timings.Add(S_timings[1]);
                        }
                        else
                        {
                            num_missed_S++;
                            S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                        }
                    }
                    else
                    {
                        num_missed_S++;
                        S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                        S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    }
                    S_timings.RemoveAt(0);
                    S_timings.RemoveAt(0);
                }
                else
                {
                    num_missed_S++;
                    S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                    S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    for (int k = 0; k < num_peaks; k++)
                    {
                        S_timings.RemoveAt(0);
                    }
                }
                expected_S1_R_diff = S1_timings[i] - R_timings[i];
                expected_S2_R_diff = S2_timings[i] - R_timings[i];
                expected_T_R_diff = T_timings_out[i] - R_timings[i];


            }

            int correct_order_errors = 0;
            int order_errors = 0;
            int final_order_errors = 0;

            if (S1_timings[0] < R_timings[0])
            {
                order_errors++;
            }
            if (T_timings_out[0] < S1_timings[0])
            {
                order_errors++;
            }
            if (S2_timings[0]  < T_timings_out[0])
            {
                order_errors++;
            }

            for (int i = 1; i < R_timings.Count()-1; i++)
            {
                if (R_timings[i] < S2_timings[i - 1])
                {
                    order_errors++;
                    R_timings[i] = S2_timings[i - 1] + (S1_timings[i] - S2_timings[i - 1]) / 2;
                    double R_temp = R_timings[i];
                    double S2_temp = S2_timings[i - 1];
                    double S1_temp = S1_timings[i];
                    if (R_timings[i] < S2_timings[i - 1])
                    {
                        correct_order_errors++;
                        R_timings[i] = S2_timings[i - 1] + 0.05;
                        if (R_timings[i] < S2_timings[i - 1])
                        {
                            final_order_errors++;
                        }
                    }
                        
                }

                if (S1_timings[i] < R_timings[i])
                {
                    S1_timings[i] = R_timings[i] + (T_timings_out[i] - R_timings[i]) / 2.0;
                    if (S1_timings[i] < R_timings[i])
                    {

                        double S1_temp = S1_timings[i];
                        double R_temp = R_timings[i];
                        double T_temp = T_timings_out[i];
                        correct_order_errors++;
                        S1_timings[i] = R_timings[i] + 0.05;
                        if (S1_timings[i] < R_timings[i])
                        {
                            final_order_errors++;
                        }
                            
                    }
                        
                    order_errors++;
                    
                }
                if (T_timings_out[i] < S1_timings[i])
                {
                    T_timings_out[i] = S1_timings[i] + (S2_timings[i] - S1_timings[i]) / 2.0;
                    if (T_timings_out[i] < S1_timings[i])
                    {
                        correct_order_errors++;
                        T_timings_out[i] = S1_timings[i] + 0.05;
                        if (T_timings_out[i] < S1_timings[i])
                        {
                            
                            final_order_errors++;
                        }
                            
                    }
                    order_errors++;

                }
                if (S2_timings[i] < T_timings_out[i])
                {
                    S2_timings[i] = T_timings_out[i] + (R_timings[i+1] - T_timings_out[i]) / 2.0;
                    if (S2_timings[i] < T_timings_out[i])
                    {
                        correct_order_errors++;
                        S2_timings[i] = T_timings_out[i] + 0.05;
                        if (S2_timings[i] < T_timings_out[i])
                        {
                            
                            final_order_errors++;
                        }
                            
                    }
                        order_errors++;
                }

                

            }

            R_timings.RemoveAt(R_timings.Count() - 1);
            if (figures)
            {
                Debug.WriteLine("Total number R timings: " + Convert.ToString(R_timings.Count()));
                Debug.WriteLine("num S estimated: " + Convert.ToString(num_missed_S));
                Debug.WriteLine("num T estimated: " + Convert.ToString(num_missed_T));
                Debug.WriteLine("num order errors: " + Convert.ToString(order_errors));
                Debug.WriteLine("num order errors after correction: " + Convert.ToString(correct_order_errors));
                Debug.WriteLine("num order errors after final correction: " + Convert.ToString(final_order_errors));
            }

             return (R_timings, T_timings_out, S1_timings, S2_timings);
        }



        static public void plot_all_timings(List<double> R_timings, List<double> T_timings, List<double> S1_timings, List<double> S2_timings, List<double> ECG, List<double> PCG, int Fs)
        {

            List<int> R_index = new List<int>();
            List<int> T_index = new List<int>();
            List<int> S1_index = new List<int>();
            List<int> S2_index = new List<int>();

            List<double> R_index_d = new List<double>();
            List<double> T_index_d = new List<double>();
            List<double> S1_index_d = new List<double>();
            List<double> S2_index_d = new List<double>();

            List<double> R_magnitude = new List<double>();
            List<double> T_magnitude = new List<double>();
            List<double> S1_magnitude = new List<double>();
            List<double> S2_magnitude = new List<double>();

            double Fs_double = Convert.ToDouble(Fs);

            for (int i = 0; i < R_timings.Count() - 1; i++)
            {
                R_index.Add(Convert.ToInt32(R_timings[i] * Fs_double));
                T_index.Add(Convert.ToInt32(T_timings[i] * Fs_double));
                S1_index.Add(Convert.ToInt32(S1_timings[i] * Fs_double));
                S2_index.Add(Convert.ToInt32(S2_timings[i] * Fs_double));

                R_magnitude.Add(ECG[R_index[i]]);
                T_magnitude.Add(ECG[T_index[i]]);
                S1_magnitude.Add(PCG[S1_index[i]]);
                S2_magnitude.Add(PCG[S2_index[i]]);

                R_index_d.Add((R_timings[i]));
                T_index_d.Add((T_timings[i]));
                S1_index_d.Add((S1_timings[i]));
                S2_index_d.Add((S2_timings[i]));
            }


            var plt = new ScottPlot.Plot(800*2, 600);
            plt.AddSignal(ECG.ToArray(), 8000, label: "ECG signal");
            plt.XLabel("time (s)");
            plt.YLabel("ECG amplitude");
            plt.Margins(0);


            var mySignalPlot2 = plt.PlotScatter(R_index_d.ToArray(), R_magnitude.ToArray(), label: "R peaks");
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            var mySignalPlot3 = plt.PlotScatter(T_index_d.ToArray(), T_magnitude.ToArray(), label: "T peaks");
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;
            plt.Legend();
            plt.Title("ECG with timings");
            plt.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\ECG_with_timings.png");


            var plt2 = new ScottPlot.Plot(800*2, 600);
            plt2.AddSignal(PCG.ToArray(), 8000, label: "PCG signal");
            plt2.YLabel("PCG amplitude");
            plt2.XLabel("time (s)");

            plt2.Margins(0);
            plt2.Title("PCG with timings");


            var mySignalPlot4 = plt2.PlotScatter(S1_index_d.ToArray(), S1_magnitude.ToArray(), label: "S1 peaks");
            mySignalPlot3.YAxisIndex = 0;
            mySignalPlot3.XAxisIndex = 0;

            var mySignalPlot5 = plt2.PlotScatter(S2_index_d.ToArray(), S2_magnitude.ToArray(), label: "S2 peaks");
            mySignalPlot4.YAxisIndex = 0;
            mySignalPlot4.XAxisIndex = 0;
            plt2.Legend();

            plt2.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\PCG_with_timings.png");
        }

        static public List<double> PanTompkinlowpass(List<double> input)
        {

            // H(z) = Y(z)/X(z) = (1 - z^-6)^2/(1-z^-1)^2
            // Y(z) (1 - z^-1)^2 = X(z) (1 - z^-6)^2
            // Y(z) (z^-2 -2z^-1 + 1) = X(z) (z^-12 - 2z^-6 + 1)
            // y[n-2] - 2y[n-1] + y[n] = x[n-12] - 2x[n-6] + x[n]
            // y[n] = 2y[n-1] - y[n-2] + x[n-12] - 2x[n-6] + x[n]

            List<double> output = new List<double>();

            for (int i = 0; i < input.Count(); i++)
            {
                if (i < 1)
                {
                    output.Add(input[i]);
                }
                else if (i < 2)
                {
                    output.Add(input[i] + 2.0 * output[i - 1]);
                }
                else if (i < 6)
                {
                    output.Add(input[i] + 2.0 * output[i - 1] - output[i - 2]);
                }
                else if (i < 12)
                {
                    output.Add(input[i] + 2.0 * output[i - 1] - output[i - 2] - 2.0 * input[i - 6]);
                }
                else
                {
                    output.Add(input[i] + 2.0 * output[i - 1] - output[i - 2] - 2.0 * input[i - 6] + input[i - 12]);
                }
                /*
                output.Add(input[i]); // + x[n]
                if (i >= 1)
                {
                    output[i] = output[i] +  2.0 * output[i - 1]; // + 2y[n-1]
                }
                if (i >= 2)
                {
                    output[i] = output[i] - output[i - 2]; // - y[n-2]
                }
                if (i >= 6)
                {
                    output[i] = output[i] - 2.0 * input[i - 6]; // -x[n-6]
                }
                if (i >= 12)
                {
                    output[i] = output[i] + input[i - 12]; // + x[n-12]
                }
                */
            }

            return output;
        }

        static public List<double> PanTompkinhighpass(List<double> input)
        {

            //  y[n] = 32x[n - 16] - y[n - 1] - x[n] + x[n - 32]

            List<double> output = new List<double>();


            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(-1.0 * input[i]);
                if (i >= 1)
                {
                    output[i] -= output[i - 1];
                }
                if (i >= 16)
                {
                    output[i] += 32.0 * input[i - 16];
                }
                if (i >= 32)
                {
                    output[i] += input[i - 32];
                }
            }

            return output;
        }

        static public List<double> PanTompkinderivative(List<double> input)
        {
            List<double> output = new List<double>();

            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i]);
                if (i >= 1)
                {
                    output[i] -= 2.0 * input[i - 1];
                }
                if (i >= 2)
                {
                    output[i] -= input[i - 2];
                }
                if (i >= 2 & i <= input.Count() - 2)
                {
                    output[i] += 2.0 * input[i + 1];
                }
                if (i >= 2 & i <= input.Count() - 3)
                {
                    output[i] += input[i + 2];
                }

                output[i] = output[i] / 8.0;
            }


            return output;
        }

        static public List<double> PanTompkinsquare(List<double> input)
        {
            List<double> output = new List<double>();

            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i] * input[i]);
            }

            return output;
        }

        static public List<double> PanTompkinint(List<double> input, int Fs)
        {
            List<double> output = new List<double>();

            int window = Convert.ToInt32(0.15 * Convert.ToDouble(Fs));
            double win_double = Convert.ToDouble(window);
            double sum = 0;


            for (int i = 0; i < window; i++)
            {
                sum += input[i] / win_double;
                output.Add(sum);
            }

            for (int i = window; i < output.Count(); i++)
            {
                sum += input[i] / win_double;
                sum -= input[i - window] / win_double;
                output.Add(sum);
            }

            return output;

        }

        static  (List<double>, List<double>) removeT(List<double> input, List<double> inputY, int Fs)
        {
            List<double> output = new List<double>();
            List<double> outputY = new List<double>();
            double min_time = 0.3;
            int min_samples = Convert.ToInt32(Convert.ToDouble(Fs) * min_time);


            for (int i = 1; i < input.Count(); i++)
            {
                if ((input[i] - input[i-1]) < min_samples)
                {
                    input.RemoveAt(i);
                    inputY.RemoveAt(i);
                }
            }
            return (input, inputY);

        }

        static List<double> findR(List<double> ECG, List<double> input, int Fs)
        {
            List<double> output = new List<double>();
            double window = 0.05;

            int window_samples = Convert.ToInt32(window * Convert.ToDouble(Fs));

            for (int i = 0; i < input.Count(); i++)
            {
                double max_value = double.MinValue;
                int max_index = 0;
                
                if (input[i] - window_samples < 0)
                {
                    for (int j = 0; j < Convert.ToInt32(input[i]) + window_samples; j++)
                    {
                        if (ECG[j] > max_value)
                        {
                            max_value = ECG[j];
                            max_index = j;
                        }

                    }
                }
                else if (input[i] + window_samples >= ECG.Count())
                {
                    for (int j = Convert.ToInt32(input[i]) - window_samples; j < ECG.Count(); j++)
                    {
                        if (ECG[j] > max_value)
                        {
                            max_value = ECG[j];
                            max_index = j;
                        }
                    }
                }
                else
                {
                    for (int j = Convert.ToInt32(input[i]) - window_samples; j < Convert.ToInt32(input[i]) + window_samples; j++)
                    {
                        if (ECG[j] > max_value)
                        {
                            max_value = ECG[j];
                            max_index = j;
                        }
                    }

                }

                
                output.Add(max_index);
            }


            return output;
        }

        static public List<double> PanTompkin(List<double> ECG, int Fs, bool figures)
        {
            // 1. Apply bandpass filter between 5 and 15 hz
            int list_length = Convert.ToInt32(Convert.ToDouble(ECG.Count()) * 200.0 / Convert.ToDouble(Fs));
            List<double> ECG_200 = Viterbi_Springer.Resample(ECG, list_length);
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG_200, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_200_resampled.png");

            }
            //List<double> ECG_low = PanTompkinlowpass(ECG_200);
            //Viterbi_Springer.PlotData(ECG_low, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_lowpass.png");
            //List<double> ECG_band = PanTompkinhighpass(ECG_low);
            //Viterbi_Springer.PlotData(ECG_band, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_highpass.png");

            List<double> ECG_band = new List<double>();
            LowpassFilterButterworthImplementation filter_Low = new LowpassFilterButterworthImplementation(15, 2, 200);
            HighpassFilterButterworthImplementation filter_high = new HighpassFilterButterworthImplementation(5, 2, 200);

            for (int i = 0; i < ECG_200.Count(); i++)
            {
                ECG_band.Add(filter_high.compute((filter_Low.compute(ECG_200[i]))));
            }
            if (figures) {
                Viterbi_Springer.PlotData(ECG_band, 200, "ECG amplitude", "time (s)", "ECG signal with 15-60 Hz band-pass filter", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_bandpass.png");

            }


            // Testing low pass
            List<double> ECG_low = PanTompkinlowpass(ECG_200);
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG_low, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_lowpass.png");
            }

            List<double> ECG_der = PanTompkinderivative(ECG_band);
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG_der, 200, "ECG amplitude", "time (s)", "ECG signal with derivative filter", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_derivative.png");

            }

            // Shift back to the left
            for (int i = 0; i < 11; i++)
            {
                ECG_der.RemoveAt(0);
                ECG_der.Add(0);
            }

            List<double> ECG_square = PanTompkinsquare(ECG_der);
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG_square, 200, "ECG amplitude", "time (s)", "filtered and squared ECG signal ", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_squared.png");
            }


            double window_duration = 0.15;
            int window = Convert.ToInt32(window_duration * 200.0);
            // ensure window is odd
            if (window % 2 == 0)
            {
                window += 1;
            }

            List<double> ECG_int = PCG_methods.smooth_data(ECG_square, window);
            PCG_methods.normalise_signal(ECG_int);
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG_int, 200, "ECG amplitude", "time (s)", "filtered, squared and normalised ECG signal with moving average filter", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_integrated.png");
            }


            if (figures)
            {
                PCG_methods.normalise_signal(ECG_200);
                PCG_methods.normalise_signal(ECG_der);

                Viterbi_Springer.MultiPlot(ECG_200, ECG_der, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\band_pass_compare.png");

            }

            // Apply thresholding. Simple thesholding should be fine 


            double threshold = 0.5;
            (List<double> R_timings, List<double> R_val) = PCG_methods.Noisy_Peak_Finding(ECG_int, 0, threshold);
            if (figures)
            {
                Viterbi_Springer.MultiPlot(ECG_int, R_timings, R_val, 200, "Squared and normalised ECG amplitude", "seconds (s)", "QRS peaks detected with thresholding", "squared, normalised and filtered ECG", "QRS peaks", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_QRS_Timings_POWER.png");
                List<double> ECG_R_vals = new List<double>();
                for (int i = 0; i < R_timings.Count(); i++)
                {
                    ECG_R_vals.Add(ECG_200[Convert.ToInt32(R_timings[i])]);


                }
                Viterbi_Springer.MultiPlot(ECG_200, R_timings, ECG_R_vals, 200, "ECG amplitude", "seconds (s)", "Identified QRS peaks on ECG signal", "ECG signal", "QRS peaks", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_QRS_Timings.png");

            }

            (R_timings, R_val) = removeT(R_timings, R_val, 200);

            if (figures)
            {
                Viterbi_Springer.MultiPlot(ECG_int, R_timings, R_val, 200, "ECG amplitude", "seconds (s)", "Identified QRS peaks on ECG signal", "squared, normalised and filtered ECG", "QRS peaks",@"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_QRS_Timings_POWER_T_Removed.png");
                List<double> ECG_R_vals = new List<double>();
                for (int i = 0; i < R_timings.Count(); i++)
                {
                    ECG_R_vals.Add(ECG_200[Convert.ToInt32(R_timings[i])]);


                }
                Viterbi_Springer.MultiPlot(ECG_200, R_timings, ECG_R_vals, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\PT_QRS_Timings_T_removed.png");

            }


            // serach for true R peak rather than just QRS. This is the highest point in the cleaned ECG sample

            R_timings = findR(ECG_200, R_timings, 200);


            return R_timings;


            //
        }

        // assuming input fs is a factor of output fs
        static public List<double> change_timing_freq(List<double> input, int input_Fs, int output_Fs)
        {
            List<double> output = new List<double>();

            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i] * (output_Fs / input_Fs));
            }

            return output;
        }

        static public List<double> shift_T_timings(List<double> input)
        {
            List<double> output = new List<double>();

            int shift_200 = 150;

            for (int i = 0; i < input.Count(); i++)
            {
                output.Add(input[i] + shift_200);
            }


            return output;
        }
    }


    



}
