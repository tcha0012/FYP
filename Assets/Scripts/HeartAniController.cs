using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using PCG_segmentation;


public class HeartAniController : MonoBehaviour
{

    
    public string FinalPath;

    public GameEngine gameEngine;
    // references to each animator
    public Animator heartCrossAni;
    public Animator bicuspidAni;
    public Animator tricuspidAni;
    public Animator rightSemiAni;
    public Animator leftSemiAni;

    // heartbeat sounds
    public AudioSource firstBeatAudio;
    public AudioSource secondBeatAudio;

    // container for the list of timings
    private TimingsContainer ecgTimings = new TimingsContainer();
    // time at which animation starts
    private float startTime;

    // variables for coroutine that need to be held between runs
    // position in the sync cycle
    private int ecgSyncIteration = 0;
    // the timing for the previous iteration
    private float ecgLastTiming = 0;

    public List<string> syncTimings = new List<string>();

    // starts the actual animation
    public void StartAnimation()
    {

        // We need to call a function to get the csv file
        /*
        Debug.Log(FileBrowser.WaitForLoadDialog(FileBrowser.PickMode.FilesAndFolders, true, null, null, "Load Files and Folders", "Load"));
        if (FileBrowser.Success) 
        {
            Debug.Log("file browser was successful");
        }
        else
        {
            Debug.Log("file browser was unsuccessful");
        }
        */
        LoadFile();
        string file_location = FinalPath;

        (List<double> R_timings, List<double> T_timings, List<double> S1_timings, List<double> S2_timings) = Segment_data.segment_data(file_location, 100, 8000);
        
        //Debug.Log(Convert.ToString(R_timings.Count());
        
        // reset all related attributes
        //ecgTimings = new TimingsContainer();
        ecgSyncIteration = 0;
        ecgLastTiming = 0;
        // populates each timings container object
        syncTimings.Clear();
        ReadDataList(syncTimings, R_timings, T_timings, S1_timings, S2_timings);
        ReadData("all_timings", ecgTimings);
        startTime = Time.time;
        StartCoroutine(EcgAnimationSync());
    }

    private void ReadDataList(List<string> syncTimings, List<double> R_timings, List<double> T_timings, List<double> S1_timings, List<double> S2_timings)
    {
        TimingsContainer timings = new TimingsContainer();
        for (int i = 0; i < R_timings.Count; i++)
        {
            
            syncTimings.Add(Convert.ToString(R_timings[i]));
            syncTimings.Add(Convert.ToString(S1_timings[i]));
            syncTimings.Add(Convert.ToString(T_timings[i]));
            syncTimings.Add(Convert.ToString(S2_timings[i]));
            
        }

        
    }

    // function that takes csv file and populates the timings container object
    private void ReadData(string dataPath, TimingsContainer timings)
    {
        // parsing csv into text asset
        TextAsset parsedData = new TextAsset();
        parsedData = Resources.Load<TextAsset>(dataPath);
        // calculates the number of columns in the inputted csv
        int columnCount = parsedData.text.Split(new string[] { "\n" }, StringSplitOptions.None)[0].Split(',').Length;
        // creates a list of strings where each item is a cell's contents
        string[] dataList = parsedData.text.Split(new string[] { ",", "\n" }, StringSplitOptions.None);
        // meta data about the csv
        int rowCount = dataList.Length / columnCount - 1;
        // for loop that populates the individual timing lists
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                timings.timingsList.Add(dataList[columnCount * (i + 1) + j]);
            }
        }
    }

    IEnumerator EcgAnimationSync()
    {
        // variables for synchronising animation
        int syncFrames = 0;
        float syncSpeed;
        //List<string> syncTimings = ecgTimings.timingsList;
        
        float timingDifference;

        // sets the sync frame and list of timings according to where in the sync cycle we are
        switch (ecgSyncIteration % 4)
        {
            case 0:
                // r at frame 12
                syncFrames = 12;
                
                break;
            case 1:
                // s1 at frame 26
                syncFrames = 14;
                firstBeatAudio.Play();
                break;
            case 2:
                // t at frame 45
                syncFrames = 19;
                break;
            case 3:
                // s2 at frame 0
                syncFrames = 5;
                secondBeatAudio.Play();
                break;
        }

        // calculates speed of animation based on sync timing
        timingDifference = float.Parse(syncTimings[ecgSyncIteration]) - ecgLastTiming;
        // sync frames divided by the timing gives the fps of the the animation divided by 60fps to derive the speed
        syncSpeed = (syncFrames / timingDifference) / 60;
        // sets the speed of the cross section
        bicuspidAni.speed = syncSpeed;
        tricuspidAni.speed = syncSpeed;
        rightSemiAni.speed = syncSpeed;
        leftSemiAni.speed = syncSpeed;
        heartCrossAni.speed = syncSpeed;
        // saves the values of current time for the next iteration
        ecgLastTiming = float.Parse(syncTimings[ecgSyncIteration]);

        // if statement ensures the animation only starts playing after the first sync speed is calculated
        if (!heartCrossAni.enabled)
        {
            bicuspidAni.enabled = true;
            tricuspidAni.enabled = true;
            rightSemiAni.enabled = true;
            leftSemiAni.enabled = true;
            heartCrossAni.enabled = true;
        }

        // increment iteration to track where in pqrst cycle we are
        ecgSyncIteration++;

        // Debug.Log(Time.time + " " + syncTimings[ecgSyncIteration]);
        // waits until the we reach the timing specified for the current iteration before running again
        float parsedTiming = float.Parse(syncTimings[ecgSyncIteration]);
        // accounts for race condition
        float adjustedTime = Time.time - startTime;
        if (adjustedTime > parsedTiming)
        {
            yield return new WaitForSeconds(timingDifference - (adjustedTime - parsedTiming));
        }
        else
        {
            yield return new WaitForSeconds(timingDifference);
        }

        // exit condition for the end of the data
        if ((ecgSyncIteration < syncTimings.Count - 1))
        {
            // continues
            StartCoroutine(EcgAnimationSync());
            yield return null;
        }
        else
        {
            // ends animation
            bicuspidAni.enabled = false;
            tricuspidAni.enabled = false;
            rightSemiAni.enabled = false;
            leftSemiAni.enabled = false;
            heartCrossAni.enabled = false;
            gameEngine.EndAni();
            yield return null;
        }
    }

    public void LoadFile()
    {
        string FileType = NativeFilePicker.ConvertExtensionToFileType("*");
        NativeFilePicker.Permission permission = NativeFilePicker.PickFile((path) =>
        {
            if (path == null)
                Debug.Log("Operation cancelled");
            else
            {
                FinalPath = path;
                Debug.Log("Picked file: " + FinalPath);
            }

        }, new string[] { FileType });
    }
}