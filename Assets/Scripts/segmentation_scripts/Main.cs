using System;
using System.Collections.Generic;
using System.Linq;
using ScottPlot;

/*

Steps to PCG segmentation

1. Take PCG data from file
2. Perform Preprocessing (bandpass filter and spike removal)
3. Get Hilbert Envelope
4. Perform Autocorrection
5. Get initial HR (using autocorrection)
6. Get Systolic time Interval 
7. Perform Segmentation 

*/
namespace PCG_segmentation
{
    class Segment_data
    {
        public static (List<double>, List<double>, List<double>, List<double>) segment_data(string data_location, double sample_length, int source_Fs)
        {
            // Creating lists to hold incoming csv data
            bool figures = false;
            bool save_output = false;


            List<int> indexes = new List<int>();
            List<double> timings = new List<double>();
            List<double> ECG = new List<double>();
            List<double> PCG = new List<double>();

            // Adding csv data to lists
            //int source_Fs = 8000;
            // sample_length has units seconds
            //double sample_length = 15;
            PCG_methods.csv_to_list(data_location, indexes, timings, ECG, PCG, sample_length, source_Fs);

            if (figures)
            {
                Viterbi_Springer.PlotData(PCG, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\Raw_PCG.png");
            }

            // getting pre-processing

            int new_Fs = 1000;
            List<double> processed_data = new List<double>();
            PCG_methods.get_hr_preprocessing(processed_data, PCG, source_Fs, new_Fs);

            (List<double> psd, int cut_data) = Viterbi_Springer.get_PSD_feature_Springer_HMM(processed_data, new_Fs, 40, 60);
            if (figures)
            {
                Viterbi_Springer.PlotData(psd, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\psd.png");

            }
            PCG_methods.normalise_signal(psd);
            if (figures)
            {
                Viterbi_Springer.PlotData(psd, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\normalised_psd.png");

            }

            int down_sampled_Fs = 150;
            int down_sampled_length = Convert.ToInt32(Convert.ToDouble(down_sampled_Fs) / Convert.ToDouble(new_Fs) * Convert.ToDouble(processed_data.Count()));
            List<double> psd_down_sampled = Viterbi_Springer.Resample(psd, down_sampled_length);

            // Test removing data 

            //List<double> cropped_processed_PCG = PCG_methods.crop_data(processed_data, new_Fs, cut_data, new_Fs);

            List<double> PCG_down_sampled = Viterbi_Springer.Resample(processed_data, down_sampled_length);

            (List<double> peak_indices, List<double> peak_values) =  PCG_methods.Noisy_Peak_Finding(psd_down_sampled, 5, 0);
            if (figures)
            {
                Viterbi_Springer.MultiPlot(psd_down_sampled, peak_indices, peak_values, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\Peaks.png");
            }


            (List<double> S_indices, List<double> S_Amplitudes) = PCG_methods.peak_post_processing(peak_indices, peak_values, down_sampled_Fs);

            if (figures)
            {
                Viterbi_Springer.MultiPlot(psd_down_sampled, S_indices, S_Amplitudes, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\Peaks_after_post_processing.png");
            }
            //Viterbi_Springer.TriplePlot(psd_down_sampled, S1_indexes, S1_Amplitudes, S2_indexes, S2_Amplitudes, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\S1_S2_segmentation.png");




            PCG_methods.normalise_signal(PCG_down_sampled);
            PCG_methods.normalise_signal(psd_down_sampled);

            if (figures)
            {
                Viterbi_Springer.MultiPlot(psd_down_sampled, PCG_down_sampled, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\psd_pcg_compare.png");
            }


            // ECG Segmentation
            if (figures)
            {
                Viterbi_Springer.PlotData(ECG, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\Raw_ECG.png");
            }


            // Cleaning data
            List<double> cleaned_ECG = ECG_methods.ecg_clean(ECG, source_Fs);
            PCG_methods.normalise_signal(cleaned_ECG);

            if (figures)
            {
                Viterbi_Springer.PlotData(cleaned_ECG, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\cleaned_ECG.png");
            }

            // I should find the R peaks before smoothing the data. As they are more pronounced at this stage.
            double R_threshhold = 2.0;
            (List<double> R_indices, List<double> R_values) = ECG_methods.ECG_peak_detection(cleaned_ECG, R_threshhold);

            if (figures)
            {
                Viterbi_Springer.MultiPlot(cleaned_ECG, R_indices, R_values, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\R_Peaks.png");
            }
            // Smoothing signal with moving average
            double window_duration = 0.1;
            int window = Convert.ToInt32(window_duration * Convert.ToDouble(source_Fs));
            // ensure window is odd
            if (window % 2 == 0)
            {
                window += 1;
            }
            List<double> smoothed_signal = PCG_methods.smooth_data(cleaned_ECG, window);

            if (figures)
            {
                Viterbi_Springer.PlotData(smoothed_signal, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\smoothed_signal.png");
            }

            // Remove the R peak so they dont interfere with peak finding of P 
            double R_window_duration = 0.05;
            int R_window = Convert.ToInt32(window_duration * Convert.ToDouble(source_Fs));
            if (R_window % 2 == 0)
            {
                R_window += 1;
            }
            List<double> smoothed_signal_without_R = ECG_methods.Remove_R_Peaks(smoothed_signal, R_indices, R_window);

            if (figures)
            {
                Viterbi_Springer.PlotData(smoothed_signal_without_R, 1, "amplitude", @"C:\github\PCG_Segmentation\PCG_segmentation\plots\no_R.png");
            }

            // Use peak finding to identify the T wave

            (List<double> T_indices, List<double> T_values) = ECG_methods.ECG_peak_detection(smoothed_signal_without_R, 0.25);

            if (figures)
            {
                Viterbi_Springer.MultiPlot(smoothed_signal_without_R, T_indices, T_values, @"C:\github\PCG_Segmentation\PCG_segmentation\plots\T_peaks.png");
            }


            List<double> S_timings = PCG_methods.samples2Timings(S_indices, down_sampled_Fs);
            //List<double> S1_timings = PCG_methods.samples2Timings(S1_indexes, down_sampled_Fs);
            //List<double> S2_timings = PCG_methods.samples2Timings(S2_indexes, down_sampled_Fs);
            List<double> pre_T_timings = PCG_methods.samples2Timings(T_indices, source_Fs);
            List<double> pre_R_timings = PCG_methods.samples2Timings(R_indices, source_Fs);

            

            (List<double> R_timings, List<double> T_timings, List<double> S1_timings, List<double> S2_timings) = ECG_methods.post_processing(pre_R_timings, pre_T_timings, S_timings);


            if (save_output)
            {
                PCG_methods.save_to_csv(S1_timings, @"C:\github\PCG_Segmentation\S1_timings.csv");
                PCG_methods.save_to_csv(S2_timings, @"C:\github\PCG_Segmentation\S2_timings.csv");
                PCG_methods.save_to_csv(T_timings, @"C:\github\PCG_Segmentation\T_timings.csv");
                PCG_methods.save_to_csv(R_timings, @"C:\github\PCG_Segmentation\R_timings.csv");
            }

            if (figures)
            {
                ECG_methods.plot_all_timings(R_timings, T_timings, S1_timings, S2_timings, ECG, PCG, source_Fs);
            }

            return (R_timings, T_timings, S1_timings, S2_timings);

        }


    }
}
