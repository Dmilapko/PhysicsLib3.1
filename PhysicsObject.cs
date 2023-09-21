using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using MonoHelper;

namespace Physics
{
    public class Force
    {
        /// <summary>
        /// Position of object to apply the force
        /// </summary>
        public PointD pos;
        /// <summary>
        /// Angle in radians
        /// </summary>
        public double angle;
        /// <summary>
        /// Force in N(Newtons)
        /// </summary>
        public double force;

        public Force(PointD _pos, double _angle, double _force)
        {
            pos = _pos;
            angle = _angle;
            force = _force;
        }
    }

    public class PhysicsObject
    {
        //public MonentOfInertia;
        public double MonentOfInertiaKoef;
        /// <summary>
        /// Angular speed in degrees per second
        /// </summary>
        public double angularvelocity = 0;
        /// <summary>
        /// Mass of object in kilos
        /// </summary>
        virtual public double mass { get; set; }
        /// <summary>
        /// Vector of speed
        /// </summary>
        public PointD speed = new PointD(0, 0);

        public PhysicsObject()
        {

        }

        public PhysicsObject(double _MomentOfInertiaKoef, double _mass)
        {
            MonentOfInertiaKoef = _MomentOfInertiaKoef;
            mass = _mass;
        }

        /// <summary>
        /// Apply force in N(Newtons)
        /// </summary>
        /// <param name="pos">Position of object to apply the force</param>
        /// <param name="angle">Angle in radians</param>
        /// <param name="force">Force in N(Newtons)</param>
        public virtual void ApplyForce(PointD pos, double angle, double force)
        {
            //angle - Math.Atan2(pos.X, pos.Y) angle between ANGLE and tangent
            double torque = (force * Math.Sqrt(pos.X * pos.X + pos.Y * pos.Y) * Math.Sin(angle - Math.Atan2(pos.X, pos.Y)));
            angularvelocity += (torque / MonentOfInertiaKoef / mass);

            speed.X += (float)(force * Math.Sin(angle) / mass);
            speed.Y += (float)(force * Math.Cos(angle) / mass);
        }

        /// <summary>
        /// Apply force in N(Newtons)
        /// </summary>
        /// <param name="force">Force applied</param>
        public virtual void ApplyForce(Force force)
        {
            ApplyForce(force.pos, force.angle, force.force);
        }
    }
}
