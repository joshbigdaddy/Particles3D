using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class SPHv2 : MonoBehaviour
{

    //Public values that can be modified in the inspector while running the script
    public float interaction_radius = 0.2f;
    public List<Particle> allParticlesInSimulation;
    public int num_particles_simulation = 1000;
    public GameObject source_component;
    public Vector2 top_right_constraint = new Vector2(0, 0);
    public Vector2 top_left_constraint = new Vector2(10, 0);
    public Vector2 bot_right_constraint = new Vector2(0, 10);
    public Vector2 bot_left_constraint = new Vector2(10, 10);
    //Our class particle, in which we store every parameter we need to use when computing our simulation
    public class Particle
    {
        public int idPart;
        public GameObject particle;
        public Vector3 velocity;
        public Vector3 gravity;
        public int hashId = 0;
        public float mass = 1f;
        public float density;
        public Vector3 acceleration;
        public List<Particle> neighbours;
        public int collisions;

        //Constructor
        public Particle()
        {

            particle = null;
            gravity = new Vector3(0, 0, 0);
            velocity = new Vector3(0, 0, 0);
            acceleration = new Vector3(0, 0, 0);
            density = 0;
            neighbours = new List<Particle>();
            collisions = 0;
        }
    }

    //Here we are going to break all the space we want in small pieces, then we are going to make a projection of our particles in a 2D Grid so as we can calculate
    //our neighbors in 2D to make the operations quicker, once we got our potential neighbors we discard all of the ones that are not optimal
    public int getPositionInHash(Particle p)
    {

        int result = 0;

        Vector2 projection = new Vector2(p.particle.transform.position.x, p.particle.transform.position.z);

        result = (int)(Mathf.Floor(projection.x / (interaction_radius * 2)) * 10 + Mathf.Floor(projection.y / (interaction_radius * 2)));

        return result;
    }

    List<Particle> findPotentialNeighbors(Particle p, Hashtable spatialHash, int positionInHash)
    {
        int[] neighbors = new int[] { -4, -3, -2, -1, 1, 2, 3, 4 };
        List<Particle> potentialNeighbors = new List<Particle>();
        for (int i = 0; i < neighbors.Length; i++)
        {
            if (spatialHash.ContainsKey(positionInHash + neighbors[i]))
            {
                potentialNeighbors.AddRange((List<Particle>)spatialHash[positionInHash + neighbors[i]]);
            }
        }

        return potentialNeighbors;

    }
    //Here we are going to take our constraints and generate diferent initial positions taking into account 
    public Vector3 initializePositionWithConstraints()
    {
        Vector3 result;
        //top_left(x,y) ================== top_right (x,y)
        // |                                     |
        // |                                     |       
        // |                                     |
        // |                                     |
        // |                                     |
        //bot_left(x,y) ================== bot_right (x,y)
        //Note: It is thought to represent a rectangle parallel to the x axis and the y axis to make everything faster and easier to calculate

        float x = Random.Range(bot_left_constraint.x, bot_right_constraint.x);
        float z = Random.Range(bot_right_constraint.y, top_right_constraint.y);
        float y = Random.Range(0, 10);
        return result = new Vector3(x, y, z);

    }
    // Use this for initialization.
    //Here we are going to populate our scene with particles and we are going to initialize everything here
    void Start()
    {

        allParticlesInSimulation = new List<Particle>();
        for (int i = 0; i < num_particles_simulation; i++)
        {
            Vector3 initializedPositionWithConstraints = initializePositionWithConstraints();//Here we get a Vector with all the initial positions we are going to give to our particles
            Particle particle = new Particle();
            particle.particle = Instantiate(source_component, initializedPositionWithConstraints, Quaternion.identity) as GameObject;
            allParticlesInSimulation.Add(particle);
        }
        Initialize(allParticlesInSimulation.ToArray());
    }

    // Update is called once per frame
    void Update()
    {
        Integrate(0.1f,allParticlesInSimulation.ToArray());
    }


    private float _softeningLengthSquared;


    public float SofteningLength { get; set; }
    

    //With this two methods we can generate something similar to a hive of bees or flies, so its very interesting to be used for insects simulations
    public void Initialize(Particle[] parts)
    {

        SofteningLength = 1;
        _softeningLengthSquared = SofteningLength * SofteningLength;


        for (int i = 0; i < parts.Length; i++)
        {
            parts[i].acceleration = Vector3.zero;

            for (int j = (i + 1); j < parts.Length; j++)
            {
                Vector3 r = parts[j].particle.transform.position - parts[i].particle.transform.position;
                float rm = r.sqrMagnitude + _softeningLengthSquared;
                Vector3 f = r / (Mathf.Sqrt(rm) * rm);
                parts[i].acceleration += f * parts[j].mass;
                parts[j].acceleration -= f * parts[i].mass;
            }
        }
    }

    public void Integrate(float dT, Particle[] bodies)
    {
        float dT2 = dT / 2.0f;

        for (int i = 0; i < bodies.Length; i++)
        {
            bodies[i].velocity += bodies[i].acceleration * dT2;
            bodies[i].particle.transform.position += bodies[i].velocity * dT;
            bodies[i].acceleration = Vector3.zero;

            for (int j = (i + 1); j < bodies.Length; j++)
            {
                Vector3 r = bodies[j].particle.transform.position - bodies[i].particle.transform.position;
                float rm = r.sqrMagnitude + _softeningLengthSquared;
                Vector3 f = r / (Mathf.Sqrt(rm) * rm);
                bodies[i].acceleration += f * bodies[j].mass;
                bodies[j].acceleration -= f * bodies[i].mass;
            }
            bodies[i].velocity += bodies[i].acceleration * dT2;
        } 
    }
}


