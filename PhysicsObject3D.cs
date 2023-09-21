using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.Xna.Framework.Input;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using MonoHelper;
using System.Data;

namespace PhysicsLib
{
    public class Force3D
    {
        /// <summary>
        /// Position of object to apply the force
        /// </summary>
        public Vector3 pos;
        /// <summary>
        /// Force in N(Newtons)
        /// </summary>
        public Vector3 direction;
        public float magnitude;

        public Force3D(Vector3 _pos, Vector3 _direction, float _magnitude)
        {
            pos = _pos;
            direction = _direction;
            magnitude = _magnitude;
        }
    }

    internal class PhysicsObject3D
    {
        //public delegate void MOIF();
        //public MOIF MonentOfInertia;
        public double MonentOfInertiaKoef;
        /// <summary>
        /// Angular speed in degrees per second
        /// </summary>
        public Vector3 angularvelocity = Vector3.Zero;
        /// <summary>
        /// Mass of object in kilos
        /// </summary>
        virtual public double mass { get; set; }
        /// <summary>
        /// Vector of speed
        /// </summary>
        public Vector3 speed = new Vector3(0, 0, 0);

        public PhysicsObject3D() 
        {
        }

        public void ApplyForce(Force3D force)
        {
            float realforce = force.direction.Z * force.magnitude;
            float radiusXZ = (float)Math.Sqrt(force.direction.X * force.direction.X + force.direction.Z * force.direction.Z);
        }
    }
}
