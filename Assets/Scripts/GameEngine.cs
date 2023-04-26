using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class GameEngine : MonoBehaviour
{
    public HeartAniController heartController;
    public EcgLineRenderer lineController;
    public GameObject heartModel;
    public GameObject startButton;
    public TextMeshProUGUI alertText;

    public void StartAni()
    {
        alertText.enabled = false;
        lineController.StartLine();
        heartController.StartAnimation();
        startButton.SetActive(false);
    }

    public void EndAni()
    {
        alertText.enabled = true;
        startButton.SetActive(true);
    }
}
