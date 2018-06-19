using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class SPHv2 : MonoBehaviour {

    public float interaction_radius = 0.2f;

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

        result = (int) (Mathf.Floor(projection.x / (interaction_radius * 2))*10 + Mathf.Floor(projection.y / (interaction_radius * 2)));

        return result;
    }
    
    List<Particle> findPotentialNeighbors(Particle p, Hashtable spatialHash, int positionInHash)
    {
        int[] neighbors =new int[] {-4 ,-3 ,-2 ,-1 ,1 ,2 ,3 ,4 };
        List<Particle> potentialNeighbors = new List<Particle>();
        for (int i=0; i<neighbors.Length;i++)
        {
            if (spatialHash.ContainsKey(positionInHash+neighbors[i]))
            {
                potentialNeighbors.AddRange((List<Particle>) spatialHash[positionInHash + neighbors[i]]);
            }
        }

        return potentialNeighbors;
    }
    // Use this for initialization
    void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
		
	}

    
        private bool _isInitialized;
        private float _softeningLengthSquared;

       
        public float SofteningLength { get; set; }

        public void Initialize(Particle[] parts)
        {
            

            _softeningLengthSquared = SofteningLength * SofteningLength;

            for (int i = 0; i < parts.Length; i++)
                parts[i].acceleration = Vector3.zero;

            for (int i = 0; i < parts.Length; i++)
            {
                for (int j = (i + 1); j < parts.Length; j++)
                {
                    Vector3 r = parts[j].particle.transform.position - parts[i].particle.transform.position;
                    float rm = r.sqrMagnitude + _softeningLengthSquared;
                    Vector3 f = r / (Mathf.Sqrt(rm) * rm);
                    parts[i].acceleration += f * parts[j].mass;
                    parts[j].acceleration -= f * parts[i].mass;
                }
            }
            _isInitialized = true;
        }

        public void Integrate(float dT, Particle[] bodies)
        {
            float dT2 = dT / 2.0f;

            for (int i = 0; i < bodies.Length; i++)
            {
                bodies[i].velocity += bodies[i].acceleration * dT2;
                bodies[i].particle.transform.position += bodies[i].velocity * dT;
                bodies[i].acceleration = Vector3.zero;
            }

            for (int i = 0; i < bodies.Length; i++)
            {
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


