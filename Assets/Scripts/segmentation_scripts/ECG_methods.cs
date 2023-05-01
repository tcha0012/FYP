using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using DSP;
using ScottPlot;

namespace PCG_segmentation
{
    internal class ECG_methods
    {

        public static string test_function_call()
        {
            return "Function call is working";
        }
        public static List<double> ecg_clean(List<double> ecg_signal, int sampling_rate)
        {
            List<double> cleaned_data = new List<double>();

            //LowpassFilterButterworthImplementation filter_low = new LowpassFilterButterworthImplementation(400, 2, Fs_new);

            // Remove baseline wander
            HighpassFilterButterworthImplementation filter_high = new HighpassFilterButterworthImplementation(0.05, 2, sampling_rate);
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
                output.Add(median_list[width/2]);
            }

            for (int i = input.Count() - width / 2; i < input.Count(); i++)
            {
                output.Add(0);
            }

            return output;
        }

        public static (List<double>, List<double>) ECG_peak_detection (List<double> input, double threshold)
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

        public static List<double> Remove_R_Peaks (List<double> input, List<double> R_indices, int window)
        {
            List<double> output = new List<double>();
            for (int i =0; i < input.Count(); i++)
            {
                output.Add(input[i]);
            }
            for (int i = 0; i < R_indices.Count(); i++)
            {
                if (Convert.ToInt32(R_indices[i]) - window < 0)
                {
                    for (int j = 0; j < Convert.ToInt32(R_indices[i]) + window; j++)
                    {
                        output[j] = 0.0;
                    }
                }
                else if (Convert.ToInt32(R_indices[i]) + window > input.Count())
                {
                    for (int j = Convert.ToInt32(R_indices[i]) - window; j < input.Count(); j++)
                    {
                        output[j] = 0.0;
                    }
                }
                else
                {
                    for (int j = Convert.ToInt32(R_indices[i]) - window; j < Convert.ToInt32(R_indices[i]) + window; j++)
                    {
                        output[j] = 0.0;
                    }
                }



            }



            return output;

        }

        static public (List<double>, List<double>, List<double>, List<double>) post_processing(List<double> R_timings, List<double> T_timings, List<double> S_timings)
        {
            // The order we want is P_wave, R_peaks, S1, T_peaks, S2
            // For now I am not using P_wave
            List<double> S1_timings = new List<double>();
            List<double> S2_timings = new List<double>();
            List<double> T_timings_out = new List<double>();
            double min_HR = 40.0; // BPM
            double max_diff = min_HR / 60.0;
            double max_HR = 200.0; // BPM

            double tolerance = 0.05;

            

            double expected_S1_R_diff = 0.075;
            double expected_S2_R_diff = 0.33 ;
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
            for (int i = 0; i < R_timings.Count()-1; i++)
            {
                
                int num_T_peaks = 0;
                
                int T_j = 0;

                
                while (T_timings[T_j] < R_timings[i + 1])
                {
                    num_T_peaks++;
                    T_j++;
                    if (T_j >= T_timings.Count())
                    {
                        break;
                    }
                }

                if (num_T_peaks == 0)
                {
                    T_timings_out.Add(R_timings[i] + expected_T_R_diff);
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

                    }
                    T_timings.RemoveAt(0);
                }
                else
                {
                    T_timings_out.Add(R_timings[i] + expected_T_R_diff);
                    for (int k = 0; k < num_T_peaks; k++)
                    {
                        T_timings.RemoveAt(0);
                    }
                }

                int num_peaks = 0;
                int j = 0;

                while (S_timings[j] < R_timings[i+1])
                {
                    num_peaks++;
                    j++;
                    if (j >= S_timings.Count())
                    {
                        break;
                    }
                }
                if (num_peaks == 0)
                {
                    S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                    S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                }
                else if (num_peaks == 1)
                {
                    if (S_timings[0] < (R_timings[i] + expected_S1_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S1_R_diff - tolerance)
                    {
                        S1_timings.Add(S_timings[0]);
                        S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    }
                    else if (S_timings[0] < (R_timings[i] + expected_S2_R_diff + tolerance) & S_timings[0] > R_timings[i] + expected_S2_R_diff - tolerance)
                    {
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
                            S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                        }
                    }
                    else
                    {
                        S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                        S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    }
                    S_timings.RemoveAt(0);
                    S_timings.RemoveAt(0);
                }
                else
                {
                    S1_timings.Add(R_timings[i] + expected_S1_R_diff);
                    S2_timings.Add(R_timings[i] + expected_S2_R_diff);
                    for (int k = 0; k < num_peaks; k++)
                    {
                        S_timings.RemoveAt(0);
                    }
                }
                


            }

            R_timings.RemoveAt(R_timings.Count() - 1);

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

            for (int i = 0; i < R_timings.Count()-1; i++)
            {
                R_index.Add(Convert.ToInt32(R_timings[i] * Fs_double));
                T_index.Add(Convert.ToInt32(T_timings[i] * Fs_double));
                S1_index.Add(Convert.ToInt32(S1_timings[i] * Fs_double));
                S2_index.Add(Convert.ToInt32(S2_timings[i] * Fs_double));

                R_magnitude.Add(ECG[R_index[i]]);
                T_magnitude.Add(ECG[T_index[i]]);
                S1_magnitude.Add(ECG[S1_index[i]]);
                S2_magnitude.Add(ECG[S2_index[i]]);

                R_index_d.Add(Convert.ToDouble(R_index[i]));
                T_index_d.Add(Convert.ToDouble(T_index[i]));
                S1_index_d.Add(Convert.ToDouble(S1_index[i]));
                S2_index_d.Add(Convert.ToDouble(S2_index[i]));
            }


            var plt = new ScottPlot.Plot(800, 600);
            plt.AddSignal(ECG.ToArray(), 1);
            plt.YLabel("amplitude");
            plt.Margins(0);
            

            var mySignalPlot2 = plt.PlotScatter(R_index_d.ToArray(), R_magnitude.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            var mySignalPlot3 = plt.PlotScatter(T_index_d.ToArray(), T_magnitude.ToArray());
            mySignalPlot2.YAxisIndex = 0;
            mySignalPlot2.XAxisIndex = 0;

            plt.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\ECG_with_timings.png");


            var plt2 = new ScottPlot.Plot(800, 600);
            plt2.AddSignal(PCG.ToArray(), 1);
            plt2.YLabel("amplitude");
            plt2.Margins(0);


            var mySignalPlot4 = plt2.PlotScatter(S1_index_d.ToArray(), S1_magnitude.ToArray());
            mySignalPlot3.YAxisIndex = 0;
            mySignalPlot3.XAxisIndex = 0;

            var mySignalPlot5 = plt2.PlotScatter(S2_index_d.ToArray(), S2_magnitude.ToArray());
            mySignalPlot4.YAxisIndex = 0;
            mySignalPlot4.XAxisIndex = 0;

            plt2.SaveFig(@"C:\github\PCG_Segmentation\PCG_segmentation\plots\PCG_with_timings.png");
        }

        

    }

}
