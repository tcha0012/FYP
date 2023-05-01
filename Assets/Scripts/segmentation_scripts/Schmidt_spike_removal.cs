using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PCG_segmentation
{
    class Schmidt_spike_removal
    {
        public static double GetMedian(List<double> sourceNumbers)
        {
            //Framework 2.0 version of this method. there is an easier way in F4        
            if (sourceNumbers == null || sourceNumbers.Count == 0)
                throw new System.Exception("Median of empty array not defined.");

            //make sure the list is sorted, but use a new array
            List<double> sortedPNumbers = new List<double>(sourceNumbers);
            sortedPNumbers.Sort();

            //get the median
            int size = sortedPNumbers.Count;
            int mid = size / 2;
            double median = (size % 2 != 0) ? (double)sortedPNumbers[mid] : ((double)sortedPNumbers[mid] + (double)sortedPNumbers[mid - 1]) / 2;
            return median;
        }

        public static void schmidt_spike_removal(List<double> original_signal, List<double> despiked_signal, int fs)
        {
            // Finding window size
            int windowsize = fs / 2;

            // Finding any samples outside of an integer number of windows
            int trailingsamples = original_signal.Count % windowsize;

            // Adding values to list
            List<double> max_vals = new List<double>();
            List<int> max_inds = new List<int>();
            List<List<double>> windows_list = new List<List<double>>();
            List<List<double>> windows_list_abs = new List<List<double>>();
            for (int i = 0; i < (original_signal.Count / windowsize); i++)
            {
                List<double> curr_window = new List<double>();
                List<double> curr_window_abs = new List<double>();
                for (int j = 0; j < windowsize; j++)
                {
                    curr_window.Add((original_signal[j + (windowsize * i)]));
                    curr_window_abs.Add(Math.Abs(original_signal[j + (windowsize * i)]));
                }
                windows_list.Add(curr_window);
                windows_list_abs.Add(curr_window_abs);
                // Find maximum value in each window
                double max_val = 0;
                int max_index = 0;
                for (int k = 0; k < curr_window_abs.Count; k++)
                {
                    if (curr_window_abs[k] > curr_window_abs[max_index])
                    {
                        max_index = k;
                        max_val = curr_window_abs[k];
                    }

                }
                max_vals.Add(max_val);
                max_inds.Add(max_index);
                
            }

       

            bool MAA_flag = true;
            int max_window_index = 0;
            while (MAA_flag)
            {
                double MAA_median = GetMedian(max_vals);  // Finding median of max_vals in each window
                MAA_flag = false;
                for (int i = 0; i < max_vals.Count(); i++) // Finding if there are samples greater than 3 times the median value that need to be removed
                {
                    if (max_vals[i] > (MAA_median * 3))
                    {
                        MAA_flag = true;
                    }
                }
                if (MAA_flag) // Removing largest spike if there is a spike that needs to be removed.
                {
                    // 1. Find window containing maximum value

                    for (int i = 0; i < max_vals.Count; i++)
                    {
                        if (max_vals[i] > max_vals[max_window_index])
                        {
                            max_window_index = i; // window and index within window can then be found as window = max_window_index, max_inds[max_window_index]
                        }
                    }

                    // find left side crossing
                    int left_crossing = 0;
                    int sign = Math.Sign(windows_list[max_window_index][max_inds[max_window_index]]);
                    for (int i = max_inds[max_window_index]; i >= 0; i--)
                    {
                        if (sign != Math.Sign(windows_list[max_window_index][i]))
                        {
                            left_crossing = i;
                            break;
                        }
                    }


                    // find right side crossing

                    int right_crossing = windowsize-1;
                    for (int i = max_inds[max_window_index]; i < windowsize; i++)
                    {
                        if (sign != Math.Sign(windows_list[max_window_index][i]))
                        {
                            right_crossing = i;
                            break;
                        }
                    }

                    // Set all values in spike to 0
                    for (int i = left_crossing; i <= right_crossing; i++)
                    {
                        windows_list[max_window_index][i] = 0.0001;
                        windows_list_abs[max_window_index][i] = 0.0001;
                        
                    }

                    // Recalculate MAAs

                    max_vals.Clear();
                    max_inds.Clear();

                    for(int i = 0; i < (original_signal.Count / windowsize); i++)
                    {
                        double max_val = 0;
                        int max_index = 0;
                        for (int k = 0; k < windowsize; k++)
                        {
                            if (windows_list_abs[i][k] > windows_list_abs[i][max_index])
                            {
                                max_index = k;
                                max_val = windows_list_abs[i][k];
                            }

                        }
                        max_vals.Add(max_val);
                        max_inds.Add(max_index);
                    }

                }
            }


            for (int i = 0; i < (original_signal.Count / windowsize); i++)
            {
                for (int j = 0; j < windowsize; j++)
                {
                    despiked_signal.Add(windows_list[i][j]);
                }
            }

            // add remaining values outside of window range

            for (int i = (original_signal.Count / windowsize) * windowsize; i < original_signal.Count; i++)
            {
                despiked_signal.Add(original_signal[i]);
            }
            
            


        }

    }
}
