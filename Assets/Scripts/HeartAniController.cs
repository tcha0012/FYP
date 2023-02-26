using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HeartAniController : MonoBehaviour
{
    // references to each animator
    public Animator heartCrossAni;
    // 
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

    private bool ecg_first_flag = true;
    private bool cusp_first_flag = true;
    private bool cusp_state_flag = true;

    private float start_time = 0;
    private float curr_time = 0;

    // Start is called before the first frame update
    void Start()
    {
        // populates each timings container object
        ReadData("TestEcgData", ecgTimings);
        ReadData("TestCuspData", cuspTimings);
        ReadData("TestSlData", slTimings);
        start_time = Time.time;
        
    }

    // Update is called once per frame
    void Update()
    {
        curr_time = Time.time;
        Debug.Log("Update: " + bicuspidAni.name);

        if (ecgSyncFlag)
        {
            StartCoroutine(EcgAnimationSync());
        }
        if (cuspSyncFlag)
        {
            StartCoroutine(CuspAnimationSync());
        }
        if (slSyncFlag)
        {
            StartCoroutine(SlAnimationSync());
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
        List<string> syncTimings = ecgTimings.timingsList;
        ecgSyncFlag = false;
        Debug.Log(heartCrossAni.GetCurrentAnimatorStateInfo(0).normalizedTime * 60);
        if ((float.Parse(syncTimings[ecgSyncIteration]) <= curr_time) || (ecg_first_flag)) {
            // variables for synchronising animation
            
            float syncFrames = 0;
            float syncSpeed;
            float timingDifference;

            // first increment of iteration
            if (!ecg_first_flag) {
            ecgSyncIteration++;
            }


            ecg_first_flag = false;
            float next_timing = float.Parse(syncTimings[ecgSyncIteration]);

            // sets the sync frame and list of timings according to where in the sync cycle we are
            switch (ecgSyncIteration % 3)
            {
                case 0:
                    // p at frame 12

                    // syncFrames = 22;
                    // going from frame 50 to 12
                    syncFrames = ((12.0f/60.0f - heartCrossAni.GetCurrentAnimatorStateInfo(0).normalizedTime)%1+1)%1 * 60f;

                    firstBeatAudio.Play();
                    break;
                case 1:
                    // r at frame 29
                    // syncFrames = 17;
                    // going from current frame to frame 29
                    syncFrames = ((29.0f/60.0f - heartCrossAni.GetCurrentAnimatorStateInfo(0).normalizedTime)%1+1)%1 * 60.0f;
                    break;
                case 2:
                    // t at frame 50    
                    // syncFrames = 21;
                    syncFrames = ((50.0f/60.0f - heartCrossAni.GetCurrentAnimatorStateInfo(0).normalizedTime)%1+1)%1 * 60.0f;
                    secondBeatAudio.Play();
                    break;
            }

            

            // checks for missing values from the segmentation


            // otherwise this block searches for the next valid timing and ensure the coroutine wait for that duration
            
            // setting flag for while loop


            // sets the appropraite wait time if the value is missing
            timingDifference = start_time + (next_timing) - curr_time;

            
            // set sync speed
            syncSpeed = (syncFrames / timingDifference) / 60;

            //Debug.Log((ecgSyncIteration % 3).ToString() + " " + syncFrames.ToString() + " " + timingDifference.ToString() + " " + syncSpeed);
            heartCrossAni.speed = syncSpeed;

            // if statement ensures the animation only starts playing after the first sync speed is calculated
            if (!heartCrossAni.enabled)
            {
                heartCrossAni.enabled = true;
            }

            // increment iteration to track where in pqrst cycle we are

            // exit condition for the end of the data
            if ((ecgSyncIteration < ecgTimings.timingsList.Count - 1))
            {
                // sets flag to true for next iteration
                ecgSyncFlag = true;
            }
            else
            {
                // disables animation
                heartCrossAni.enabled = false;
            }
        }

        else {
            yield return new WaitForSeconds(0.0f);
        }

        ecgSyncFlag = true;

    }

    IEnumerator CuspAnimationSync()
    {
        cuspSyncFlag = false;
        List<string> syncTimings = cuspTimings.timingsList;
        float state = bicuspidAni.GetCurrentAnimatorStateInfo(0).normalizedTime;
        
        if (state > 0.7)
        {
            cusp_state_flag = true;
        }
        //Debug.Log(bicuspidAni.GetCurrentAnimatorStateInfo(0).normalizedTime * 60);
        //if ((float.Parse(syncTimings[cuspSyncIteration]) <= curr_time) || cusp_first_flag || (cusp_state_flag && (state < 0.1))) {
        if ((cusp_state_flag && (state < 0.7)))
        {

            // variables for synchronising animation
            cusp_state_flag = false;


            float syncFrames = 60.0f;
            float syncSpeed;
            float timingDifference;

            // first increment of iteration
            if (!cusp_first_flag) {
                cuspSyncIteration++;
            }
            

            cusp_first_flag = false;
            float next_timing = float.Parse(syncTimings[cuspSyncIteration]);
            // sets the sync frame and list of timings according to where in the sync cycle we are


            // checks for missing values from the segmentation


            // otherwise this block searches for the next valid timing and ensure the coroutine wait for that duration
            
            // setting flag for while loop


            // sets the appropraite wait time if the value is missing
            timingDifference = start_time + (next_timing) - curr_time;
            // set sync speed
            
            if ((state == 0) | ((state > 0.8) & state < 1.0))
            {
                syncFrames = (( 0 -state) % 1 + 1) % 1 * 60 + 60;
            }
            else
            {
                syncFrames = ((0 -state) % 1 + 1) % 1 * 60;
            }
            syncSpeed = (syncFrames / timingDifference) / 60;

            tricuspidAni.speed = syncSpeed;
            bicuspidAni.speed = syncSpeed;

            Debug.Log(state.ToString() + " " + syncFrames.ToString() + " " + syncSpeed.ToString());

            // if statement ensures the animation only starts playing after the first sync speed is calculated
            if (!bicuspidAni.enabled && !tricuspidAni.enabled)
            {
                bicuspidAni.enabled = true;
                tricuspidAni.enabled = true;
            }

            // increment iteration to track where in pqrst cycle we are

            // exit condition for the end of the data
            if ((cuspSyncIteration < cuspTimings.timingsList.Count - 1))
            {
                // sets flag to true for next iteration
                cuspSyncFlag = true;
            }
            else
            {
                // disables animation
                bicuspidAni.enabled = false;
                tricuspidAni.enabled = false;
            }
        }

        else {
            yield return new WaitForSeconds(0.0f);
        }
        cuspSyncFlag = true;
    }

    IEnumerator SlAnimationSync()
    {
        slSyncFlag = false;
        // variables for synchronising animation
        int syncFrames = 0;
        float syncSpeed;
        List<string> syncTimings = slTimings.timingsList;
        float timingDifference;

        // sets the sync frame and list of timings according to where in the sync cycle we are
        switch (slSyncIteration % 2)
        {
            case 0:
                syncFrames = 26;
                break;
            case 1:
                syncFrames = 34;
                break;
        }

        // checks for missing values from the segmentation
        if (syncTimings[slSyncIteration] != "NULL")
        {
            // calculates speed of animation based on sync timing
            timingDifference = start_time + float.Parse(syncTimings[slSyncIteration]) - curr_time;
            // sync frames divided by the timing gives the fps of the the animation divided by 60fps to derive the speed
            syncSpeed = (syncFrames / timingDifference) / 60;
            // sets the speed of the cross section
            rightSemiAni.speed = syncSpeed;
            leftSemiAni.speed = syncSpeed;
            // saves the values of current time for the next iteration
            slLastTiming = float.Parse(syncTimings[slSyncIteration]);
        }
        else
        {
            // otherwise this block searches for the next valid timing and ensure the coroutine wait for that duration
            // first increment of iteration
            slSyncIteration++;
            // setting flag for while loop
            bool whileFlag = false;
            // initialise empty variable
            string nextValidTiming = "";
            // while loop that finds the next valid timing
            while (!whileFlag)
            {
                nextValidTiming = syncTimings[slSyncIteration];
                if (nextValidTiming != "NULL")
                {
                    whileFlag = true;
                }
                else
                {
                    slSyncIteration++;
                }
            }
            // sets the appropraite wait time if the value is missing
            timingDifference = start_time + float.Parse(nextValidTiming) - curr_time;
        }

        // if statement ensures the animation only starts playing after the first sync speed is calculated
        if (!rightSemiAni.enabled && !leftSemiAni.enabled)
        {
            rightSemiAni.enabled = true;
            leftSemiAni.enabled = true;
        }

        // increment iteration to track where in pqrst cycle we are
        slSyncIteration++;

        // waits until the we reach the timing specified for the current iteration before running again
        yield return new WaitForSeconds(timingDifference);

        // exit condition for the end of the data
        if ((slSyncIteration < slTimings.timingsList.Count - 1))
        {
            // sets flag to true for next iteration
            slSyncFlag = true;
        }
        else
        {
            // disables animation
            rightSemiAni.enabled = false;
            leftSemiAni.enabled = false;
        }
    }
}
