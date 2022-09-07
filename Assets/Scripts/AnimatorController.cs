using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AnimatorController : MonoBehaviour
{
    public Animator heartCrossAni;
    public Animator bicuspidAni;
    public Animator tricuspidAni;
    public Animator rightSemiAni;
    public Animator leftSemiAni;

    // Start is called before the first frame update
    void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {
        heartCrossAni.speed += 0.0001f;
        bicuspidAni.speed += 0.0001f;
        tricuspidAni.speed += 0.0001f;
        rightSemiAni.speed += 0.0001f;
        leftSemiAni.speed += 0.0001f;
    }
}
