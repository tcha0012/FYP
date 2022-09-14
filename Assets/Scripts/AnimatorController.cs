using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AnimatorController : MonoBehaviour
{
    // references to each animator
    public Animator heartCrossAni;
    public Animator bicuspidAni;
    public Animator tricuspidAni;
    public Animator rightSemiAni;
    public Animator leftSemiAni;

    // container for the list of timings
    private TimingsContainer timings = new TimingsContainer();

    // variables for coroutine that need to be held between runs
    // position in the sync cycle
    private int syncIteration = 0;
    // the timing for the previous iteration
    private float lastTiming = 0;
    // flag for routine wait time
    private bool syncFlag = true;

    // Start is called before the first frame update
    void Start()
    {
        // populates the timings container object
        ReadData("TestData");
    }

    // Update is called once per frame
    void Update()
    {
        if (syncFlag)
        {
            StartCoroutine(AnimationSync());
        }
    }

    // function that takes csv file and populates the timings container object
    private void ReadData(string dataPath)
    {
        // parsing csv into text asset
        TextAsset parsedData = new TextAsset();
        parsedData = Resources.Load<TextAsset>(dataPath);
        // creates a list of strings where each item is a cell's contents
        string[] dataList = parsedData.text.Split(new string[] { ",", "\n" }, StringSplitOptions.None);
        // meta data about the csv
        int columnCount = 5;
        int rowCount = dataList.Length / columnCount - 1;
        // for loop that populates the individual timing lists
        for (int i = 0; i < rowCount; i++)
        {
            timings.timingsList.Add(dataList[columnCount * (i + 1)]);
            timings.timingsList.Add(dataList[columnCount * (i + 1) + 1]);
            timings.timingsList.Add(dataList[columnCount * (i + 1) + 2]);
            timings.timingsList.Add(dataList[columnCount * (i + 1) + 3]);
            timings.timingsList.Add(dataList[columnCount * (i + 1) + 4]);
        }
    }

    IEnumerator AnimationSync()
    {
        syncFlag = false;
        // variables for synchronising animation
        int syncFrames = 12;
        float syncSpeed;
        List<string> syncTimings = timings.timingsList;
        float timingDifference;

        // sets the sync frame and list of timings according to where in the sync cycle we are
        switch (syncIteration % 5)
        {
            case 0:
                // p at frame 12
                if (syncIteration != 0)
                {
                    syncFrames = 22;
                }
                else
                {
                    syncFrames = 12;
                }
                break;
            case 1:
                // q at frame 27
                syncFrames = 15;
                break;
            case 2:
                // r at frame 29
                syncFrames = 2;
                break;
            case 3:
                // s at frame 31
                syncFrames = 2;
                break;
            case 4:
                // t at frame 50    
                syncFrames = 19;
                break;
        }

        // checks for missing values from the segmentation
        if (syncTimings[syncIteration] != "NULL")
        {
            // calculates speed of animation based on sync timing
            timingDifference = float.Parse(syncTimings[syncIteration]) - lastTiming;
            // sync frames divided by the timing gives the fps of the the animation divided by 60fps to derive the speed
            syncSpeed = (syncFrames / timingDifference) / 60;
            // sets the speed of the cross section
            heartCrossAni.speed = syncSpeed;
            // saves the values of current time for the next iteration
            lastTiming = float.Parse(syncTimings[syncIteration]);
        }
        else
        {
            // otherwise this block searches for the next valid timing and ensure the coroutine wait for that duration
            // first increment of iteration
            syncIteration++;
            // setting flag for while loop
            bool whileFlag = false;
            // initialise empty variable
            string nextValidTiming = "";
            // while loop that finds the next valid timing
            while (!whileFlag)
            {
                nextValidTiming = syncTimings[syncIteration];
                if (nextValidTiming != "NULL")
                {
                    whileFlag = true;
                }
                else
                {
                    syncIteration++;
                }
            }
            // sets the appropraite wait time if the value is missing
            timingDifference = float.Parse(nextValidTiming) - lastTiming;
        }

        // if statement ensures the animation only starts playing after the first sync speed is calculated
        if (!heartCrossAni.enabled)
        {
            heartCrossAni.enabled = true;
            // bicuspidAni.enabled = true;
            // tricuspidAni.enabled = true;
            // rightSemiAni.enabled = true;
            // leftSemiAni.enabled = true;
        }

        // increment iteration to track where in pqrst cycle we are
        syncIteration++;

        // waits until the we reach the timing specified for the current iteration before running again
        yield return new WaitForSeconds(timingDifference);

        // exit condition for the end of the data
        if ((syncIteration < timings.timingsList.Count - 1))
        {
            // sets flag to true for next iteration
            syncFlag = true;
        }
        else
        {
            // disables animation
            heartCrossAni.enabled = false;
            // bicuspidAni.enabled = false;
            // tricuspidAni.enabled = false;
            // rightSemiAni.enabled = false;
            // leftSemiAni.enabled = false;
        }
    }
}
