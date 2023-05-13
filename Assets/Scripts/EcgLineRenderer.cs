using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

[RequireComponent(typeof(RectTransform))]
public class EcgLineRenderer : MonoBehaviour
{
    [Range(.1f, 10f)]
    public float thickness = 1f;    // Thickness of the line renderer
    public Color color = Color.red; // Color of the line renderer

    public bool traceInitialLine = true;    // Trace line at start
    public float yMin = 0f;                 // Data minimum value
    public float yMax = 100f;               // Data maximum value
    public float xRange = 10f;              // Time period shown

    LineRenderer lineRenderer;  // Inner component tracing the line
    RectTransform T;            // Inner component holding the canvas transform

    private TimingsContainer lineTimings = new TimingsContainer();
    private int valueIndex = 0;

    public List<double> ECG_data = new List<double>();

    // Called when application or editor opens
    void Awake()
    {
        ReadData("PCG_500_shortened", lineTimings);
        InitInnerComponents();
        UpdateLineProperties();

        // Show line placeholder in editor only
        if (!Application.isPlaying)
            TracePlaceholderLine();
    }

    // Called when the component is being enabled
    void OnEnable()
    {
        lineRenderer.enabled = true;
    }

    // Called when the component is being disabled
    void OnDisable()
    {
        lineRenderer.enabled = false;
    }

    public void prepare(string path)
    {
        List<double> ECG_500 = getECG(path, 60, 8000, 500);

    }

    public List<double> getECG(string path, double length, int Fs, int output_Fs)
    {
        using (var reader = new StreamReader(path))
        {
            List<double> indexes = new List<double>();
            List<double> timings = new List<double>();
            List<double> ECG = new List<double>();
            List<double> PCG = new List<double>();

            int i = 0;
            bool first_row_flag = true;
            while ((!reader.EndOfStream) && (i < Convert.ToInt32(Convert.ToDouble(Fs) * length)))
            {

                var line = reader.ReadLine();
                var values = line.Split(';', ',');

                if (!first_row_flag)
                {
                    indexes.Add(Convert.ToInt32(values[0]));
                    timings.Add(Convert.ToDouble(values[1]));
                    ECG.Add(Convert.ToDouble(values[2]));
                    PCG.Add(Convert.ToDouble(values[3]));
                }
                first_row_flag = false;
                i++;
            }

            List<double> ECG_500 = new List<double>();
            for (int j = 0; j < (Fs / output_Fs); j++)
            {
                ECG_500.Add(j * (Fs / output_Fs));
            }

            return ECG_500;
        }
    }

    public void StartLine()
    {
        valueIndex = 0;
        // adjust here and line 122 for time between points
        InvokeRepeating("DrawLine", 0.0f, 0.002f);
    }

    void DrawLine()
    {

        if (valueIndex < lineTimings.timingsList.Count)
        {
            AddPoint(float.Parse(lineTimings.timingsList[valueIndex]));
            valueIndex++;
        }
        else
        {
            CancelInvoke();
        }
    }

    void DrawLine2()
    {

        if (valueIndex < ECG_data.Count)
        {
            AddPoint((ECG_data[valueIndex]));
            valueIndex++;
        }
        else
        {
            CancelInvoke();
        }
    }

    // Initializes line and transform components
    void InitInnerComponents()
    {
        // Transform
        T = (RectTransform)this.transform;

        // Line
        lineRenderer = gameObject.GetComponent<LineRenderer>();
        if (!lineRenderer)
        {
            lineRenderer = gameObject.AddComponent<LineRenderer>();
            lineRenderer.material = new Material(Shader.Find("Sprites/Default"));
        }
        lineRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
        lineRenderer.receiveShadows = false;
        lineRenderer.allowOcclusionWhenDynamic = false;
        lineRenderer.useWorldSpace = false; // use canvas relative space
        lineRenderer.alignment = LineAlignment.TransformZ; // orthogonal to canvas
        lineRenderer.sortingOrder = 3;  // in front of image canvas (?)
        lineRenderer.positionCount = 0; // no points
    }

    // Show line placeholder
    void TracePlaceholderLine()
    {
        if (!lineRenderer)
            return;

        // Remove any points
        lineRenderer.positionCount = 0;

        // Trace a flat line in the middle
        var mean = (yMax + yMin) / 2;
        AddPoint(mean);
    }

    // Update line color and width
    void UpdateLineProperties()
    {
        if (!lineRenderer)
            return;

        lineRenderer.widthMultiplier = thickness / 100 * T.rect.height;
        lineRenderer.startColor = color;
        lineRenderer.endColor = color;
    }

    // Append data to the line
    void AddPoint(double val)
    {
        // Add a point on the left
        lineRenderer.positionCount++;

        double dx = DeltaTimeToDeltaX(0.002);

        // Move previous position to the left
        for (int i = lineRenderer.positionCount - 1; i > 0; --i)
        {
            var pos = lineRenderer.GetPosition(i - 1);
            pos.x = pos.x - (float)dx;
            lineRenderer.SetPosition(i, pos);
        }

        // Put new point on far right (far right, id = 0)
        float x = (1 - T.pivot.x) * T.rect.width;
        float y = ValueToY(Mathf.Clamp((float)val, yMin, yMax));
        lineRenderer.SetPosition(0, new Vector3(x, y));

        // Check if points out of bounds need to be removed
        float originX = -T.pivot.x * T.rect.width;
        for (int i = lineRenderer.positionCount - 1; i > 0; --i)
        {
            // If point not out of bounds, break
            var pos1 = lineRenderer.GetPosition(i);
            if (pos1.x >= originX)
                break;

            // If next point in bound, move current point to the intersection
            // of the line with the left bound of the canvas
            var pos2 = lineRenderer.GetPosition(i - 1);
            if (pos2.x > originX)
            {
                // Calculate slop-intercept between those two points
                var m = (pos2.y - pos1.y) / (pos2.x - pos1.x);
                var b = pos1.y - m * pos1.x;

                // Place last point on the far left
                pos1.x = originX;
                pos1.y = m * originX + b;
                lineRenderer.SetPosition(i, pos1);
                break;
            }

            // If next point out of bounds too, decrease and go to next
            lineRenderer.positionCount--;
        }

        // Add point of far left if we only have one point on far right
        if (lineRenderer.positionCount == 1 && traceInitialLine)
        {
            lineRenderer.positionCount++;
            lineRenderer.SetPosition(1, new Vector3(originX, y));
        }
    }

    double DeltaTimeToDeltaX(double t)
    {
        double p = t / xRange;
        return p * T.rect.width;
    }


    float ValueToY(float val)
    {
        float p = (val - yMin) / (yMax - yMin);
        return (p - T.pivot.y) * T.rect.height;
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
}
