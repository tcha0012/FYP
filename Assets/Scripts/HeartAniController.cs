using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HeartAniController : MonoBehaviour
{
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
    private TimingsContainer cuspTimings = new TimingsContainer();
    private TimingsContainer slTimings = new TimingsContainer();

    // variables for coroutine that need to be held between runs
    // position in the sync cycle
    private int ecgSyncIteration = 0;
    private int cuspSyncIteration = 0;
    private int slSyncIteration = 0;
    // the timing for the previous iteration
    private float ecgLastTiming = 0;
    private float cuspLastTiming = 0;
    private float slLastTiming = 0;
    // flag for routine wait time
    private bool ecgSyncFlag = true;
    private bool cuspSyncFlag = true;
    private bool slSyncFlag = true;

    // Start is called before the first frame update
    void Start()
    {
        // populates each timings container object
        ReadData("TestEcgData", ecgTimings);
        ReadData("TestCuspData", cuspTimings);
        ReadData("TestSlData", slTimings);
        StartCoroutine(EcgAnimationSync());
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
        ecgSyncFlag = false;
        // variables for synchronising animation
        int syncFrames = 0;
        float syncSpeed;
        List<string> syncTimings = ecgTimings.timingsList;
        float timingDifference;

        // sets the sync frame and list of timings according to where in the sync cycle we are
        switch (ecgSyncIteration % 3)
        {
            case 0:
                // p at frame 12
                if (ecgSyncIteration != 0)
                {
                    syncFrames = 20;
                }
                else
                {
                    syncFrames = 10;
                }
                // firstBeatAudio.Play();
                break;
            case 1:
                syncFrames = 10;
                break;
            case 2:
                syncFrames = 10;
                // secondBeatAudio.Play();
                break;
            case 3:
                syncFrames = 10;
                // secondBeatAudio.Play();
                break;
            case 4:
                syncFrames = 10;
                // secondBeatAudio.Play();
                break;
        }

        // checks for missing values from the segmentation
        if (syncTimings[ecgSyncIteration] != "NULL")
        {
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
        }
        else
        {
            // otherwise this block searches for the next valid timing and ensure the coroutine wait for that duration
            // first increment of iteration
            ecgSyncIteration++;
            // setting flag for while loop
            bool whileFlag = false;
            // initialise empty variable
            string nextValidTiming = "";
            // while loop that finds the next valid timing
            while (!whileFlag)
            {
                nextValidTiming = syncTimings[ecgSyncIteration];
                if (nextValidTiming != "NULL")
                {
                    whileFlag = true;
                }
                else
                {
                    ecgSyncIteration++;
                }
            }
            // sets the appropraite wait time if the value is missing
            timingDifference = float.Parse(nextValidTiming) - ecgLastTiming;
        }

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

        // waits until the we reach the timing specified for the current iteration before running again
        yield return new WaitForSeconds(timingDifference);

        // exit condition for the end of the data
        if ((ecgSyncIteration < ecgTimings.timingsList.Count - 1))
        {
            // sets flag to true for next iteration
            StartCoroutine(EcgAnimationSync());
        }
        else
        {
            bicuspidAni.enabled = false;
            tricuspidAni.enabled = false;
            rightSemiAni.enabled = false;
            leftSemiAni.enabled = false;
            heartCrossAni.enabled = false;
            yield return null;
        }
    }
}