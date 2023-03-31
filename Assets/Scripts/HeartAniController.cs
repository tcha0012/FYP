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

    // variables for coroutine that need to be held between runs
    // position in the sync cycle
    private int ecgSyncIteration = 0;
    // the timing for the previous iteration
    private float ecgLastTiming = 0;

    // Start is called before the first frame update
    void Start()
    {
        // populates each timings container object
        ReadData("all_timings", ecgTimings);
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
        // variables for synchronising animation
        int syncFrames = 0;
        float syncSpeed;
        List<string> syncTimings = ecgTimings.timingsList;
        float timingDifference;

        // sets the sync frame and list of timings according to where in the sync cycle we are
        switch (ecgSyncIteration % 5)
        {
            case 0:
                // p at frame 12
                syncFrames = 12;
                firstBeatAudio.Play();
                break;
            case 1:
                // s at frame 22
                syncFrames = 10;
                break;
            case 2:
                // s1 at frame 26
                syncFrames = 4;
                break;
            case 3:
                // t at frame 45
                syncFrames = 19;
                break;
            case 4:
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

        Debug.Log(Time.time + " " + syncTimings[ecgSyncIteration]);
        // waits until the we reach the timing specified for the current iteration before running again
        float parsedTiming = float.Parse(syncTimings[ecgSyncIteration]);
        if (Time.time > parsedTiming)
        {
            yield return new WaitForSeconds(timingDifference - (Time.time - parsedTiming));
        }
        else
        {
            yield return new WaitForSeconds(timingDifference);
        }

        // exit condition for the end of the data
        if ((ecgSyncIteration < ecgTimings.timingsList.Count - 1))
        {
            StartCoroutine(EcgAnimationSync());
            yield return null;
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