using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using MonoHelper;

namespace Physics
{
    public class MatrixPhysics
    {
        private class Collision
        {
            public PointD pos;
            public double surface_angle;
            public double parallel_speed;
            public double perpendicular_speed;
            public double object_rotation;

            public Collision(PointD _pos, double _surface_angle, double _perpendicular_speed, double _parallel_speed, double _object_rotation)
            {
                pos = _pos;
                surface_angle = _surface_angle;
                parallel_speed = _parallel_speed;
                perpendicular_speed = _perpendicular_speed;
                object_rotation = _object_rotation;
            }
        }

        public class MP_Object : PhysicsObject
        {
            /// <summary>
            /// Physics points in PIXELS
            /// </summary>
            public List<PointD> hitpoints = new List<PointD>();
            /// <summary>
            /// Death points in PIXELS
            /// </summary>
            public List<PointD> deathpoints = new List<PointD>();
            public PointD position;
            public double rotation = 0;
            public List<Force> forces = new List<Force>();
            public bool alive = true;
            public event EventHandler on_death = null;
            public event EventHandler on_collision = null;
            public double friction_coefficient;
            public bool active = true;
            public List<bool> ever_collided = new List<bool>();
            public List<bool> now_collided = new List<bool>();

            public MP_Object(double _MomentOfInertiaKoef, double _mass, double _friction_coefficient, PointD _position)
            {
                MonentOfInertiaKoef = _MomentOfInertiaKoef;
                mass = _mass;
                position = _position;
                friction_coefficient = _friction_coefficient;
                InitializePoints();
                InitializeDeathPoints();
                for (int i = 0; i < hitpoints.Count; i++) ever_collided.Add(false);
                for (int i = 0; i < hitpoints.Count; i++) now_collided.Add(false);
            }

            public virtual void InitializePoints()
            {
                List<PointD> points = new List<PointD>();
            }

            public virtual void InitializeDeathPoints()
            {
                List<PointD> deathpoints = new List<PointD>();
            }

            public void ApplyForce(PointD pos, double angle, double force, bool save = true)
            {
                base.ApplyForce(pos, angle, force);
                if (save) forces.Add(new Force(pos, angle, force));
            }

            public virtual void Run()
            {

            }
            public virtual void BeforeRun()
            {

            }
            public virtual void AfterRun()
            {

            }

            public void Kill()
            {
                alive = false;
                on_death?.Invoke(null, null);
            }

            public void Collide()
            {
                on_collision?.Invoke(null, null);
            }
        }

        public List<MP_Object> objects;
        public List<double> matrix = new List<double>();
        public double pim;
        public double fps;
        public int precision;

        public void SetMatrix(List<List<bool>> _matrix)
        {
            matrix.Clear();
            for (int x = 0; x < _matrix.Count; x++)
            {
                int height;
                for (height = 0; height < _matrix[x].Count; height++)
                {
                    if (!_matrix[x][height]) break;
                }
                matrix.Add(height);
            }
        }

        public void SetMatrix(List<double> _matrix)
        {
            matrix = new List<double>(_matrix);
        }

        public void SetMatrix(Microsoft.Xna.Framework.Graphics.Texture2D texture, List<Microsoft.Xna.Framework.Color> true_data)
        {
            Microsoft.Xna.Framework.Color[] cd = new Microsoft.Xna.Framework.Color[(uint)(texture.Width * texture.Height)];
            texture.GetData(cd);
            matrix.Clear();
            for (int x = 0; x < texture.Width; x++)
            {
                int height;
                for (height = 0; height < texture.Height; height++)
                {
                    if (!cd[height * texture.Width + x].In(true_data.ToArray())) break;
                }
                matrix.Add(height);
            }
        }

        public Microsoft.Xna.Framework.Graphics.Texture2D GetTexture(Microsoft.Xna.Framework.Graphics.GraphicsDevice graphics, int height, Microsoft.Xna.Framework.Color surface_color, Microsoft.Xna.Framework.Color outline_color, bool smoothed = true)
        {
            Microsoft.Xna.Framework.Graphics.Texture2D texture = new Microsoft.Xna.Framework.Graphics.Texture2D(graphics, matrix.Count, height);
            Microsoft.Xna.Framework.Color[] cd = new Microsoft.Xna.Framework.Color[height * matrix.Count];
            texture.GetData(cd);
            for (int x = 0; x < matrix.Count; x++)
            {
                int pixel_h = (int)matrix[x];
                for (int y = 0; y < pixel_h - 1; y++)
                {
                    if (((x != 0) && (matrix[x - 1] <= y)) || ((x != matrix.Count - 1) && (matrix[x + 1] <= y))) cd[y * matrix.Count + x] = outline_color;
                    else cd[y * matrix.Count + x] = surface_color;
                }
                cd[(pixel_h - 1) * matrix.Count + x] = MHeleper.MixTwoColorsNA(new Microsoft.Xna.Framework.Color(surface_color, (float)matrix[x] - pixel_h), new Microsoft.Xna.Framework.Color(outline_color, 1 - (float)matrix[x] + pixel_h));
                cd[(pixel_h) * matrix.Count + x] = new Microsoft.Xna.Framework.Color(outline_color, (float)matrix[x] - pixel_h);
            }
            if (smoothed)
            {
                for (int x = 1; x < matrix.Count; x++)
                {
                    if (Math.Abs(matrix[x-1]-matrix[x])>=1)
                    {
                        int px = x;
                        int side_x = x - 1;
                        if (matrix[x] > matrix[x - 1]) { px--; side_x++; }
                        int top = (int)Math.Max(matrix[x], matrix[x - 1]);
                        int bottom = (int)Math.Min(matrix[x], matrix[x - 1]);
                        for (int y = bottom; y < top; y++)
                        {
                            double coef = (top - y - 0.5) / (double)(top - bottom);
                            cd[y * matrix.Count + px] = new Microsoft.Xna.Framework.Color(outline_color, (float)coef);
                            if (y!=top-1)cd[y * matrix.Count + side_x] = MHeleper.MixTwoColors(surface_color, new Microsoft.Xna.Framework.Color(outline_color, 1 -(float)coef));
                        }
                    }
                }
                for (int x = 0; x < matrix.Count; x++)
                {
                    bool found = false;
                    for (int y = height - 1; y != -1; y--)
                    {
                        if (found)
                        {
                            cd[y * matrix.Count + x] = surface_color;
                        }
                        else if (cd[y * matrix.Count + x] == outline_color) found = true;
                    }
                }

            }
            texture.SetData(cd);
            return texture;
        }

        static public List<double> SmoothSurface(List<int> list)
        {
            int prev_id = 0;
            List<double> res = new List<double>();
            list.Add(list[list.Count - 1] + 1);
            for (int i = 1; i < list.Count; i++)
            {
                if (list[i] != list[prev_id])
                {
                    for (int j = prev_id; j < i; j++)
                    {
                        res.Add(list[prev_id] + (j - prev_id) / (double)(i - prev_id) * (list[i] - list[prev_id]));
                    }
                    prev_id = i;
                }
            }
            return res;
        }

        public void SmothSurface()
        {
            SetMatrix(SmoothSurface(matrix.Select(i => (int)i).ToList()));
        }

        static public List<double> CreateSurface(int width, int min_width_diapason, int max_width_diapason, int min_height, int max_height)
        {
            List<Microsoft.Xna.Framework.Vector2> points = new List<Microsoft.Xna.Framework.Vector2>();
            double prev = min_height + MHeleper.RandomDouble() * (max_height - min_height);
            int x = 0, rd =  0;
            //for (int x = 0; x - rd < width; x += rd)
            do {
                double totop = max_height - prev;
                double tobottom = prev - min_height;
                double koef = 0;
                if (totop < 100) koef = (-1 + (totop / 100f)) / 2;
                if (tobottom < 100) koef = (1 - (tobottom / 100f)) / 2;
                double height = (double)(prev + (MHeleper.RandomDouble() - 0.5 + koef) * 200);
                if (height < min_height) min_height = 100;
                if (height > max_height) height = max_height;
                prev = height;
                points.Add(new Microsoft.Xna.Framework.Vector2(x, (float)height));
                rd = (int)Math.Round(min_width_diapason + MHeleper.RandomDouble() * (max_width_diapason - min_width_diapason));
                x += rd;
            } while (x - rd < width);
            return SmoothSurface(MHeleper.CreateCurve(width, max_height, points));
        }

        public MatrixPhysics(List<MP_Object> _objects, double pixels_in_meter, double _fps, int _precision)
        {
            objects = _objects;
            pim = pixels_in_meter;
            fps = _fps;
            precision = _precision;
        }

        public PointD GetMatrixPosition(PointD point, PointD position, double rotation)
        {
            double l = point.Length();
            double d = rotation + (double)Math.Atan2(point.X, point.Y);
            return new PointD(position.X + Math.Sin(d) * l, position.Y + Math.Cos(d) * l);
        }

        public bool? GetMatrixState(PointD position)
        {
            position *= pim;
            if ((position.X > -0.5) && (position.X < (matrix.Count - 0.5))) 
            {
                return position.Y < GetMatrixHeight(position.X);
            }
            else return null;
        }

        public double GetMatrixHeight(double x)
        {
            int prevpos = Math.Max(0, (int)x);
            int nextpos = Math.Min(matrix.Count - 1, (int)x + 1);
            return matrix[prevpos] + (x - prevpos) / 1d * (matrix[nextpos] - matrix[prevpos]);
        }

    /*    public double GetSurfaceAngle(int position)
        {
            int toend = matrix.Count() - position - 1;
            PointD point;
            if (position < 10) point = new PointD(20, matrix[20] - matrix[0]);
            else if (toend < 10 ) point = new PointD(20, matrix[matrix.Count - 1] - matrix[matrix.Count - 21]);
            else point = new PointD(20, matrix[position + 10] - matrix[position - 10]);
            return Math.Atan2(point.X, point.Y);

        }*/

        public double GetSurfaceAngle(double position)
        {
            double toend = matrix.Count() - position - 1;
            PointD point;
            if (position < 10) point = new PointD(21, matrix[20] - matrix[0]);
            else if (toend < 10) point = new PointD(21, matrix[matrix.Count - 1] - matrix[matrix.Count - 21]);
            else point = new PointD(21, GetMatrixHeight(position + 10) - GetMatrixHeight(position - 10));
            return Math.Atan2(point.X, point.Y);

        }

        /*private void RunObjPhysics(MP_Object obj)
        {
            if (obj.alive)
            {
                obj.forces.Clear();
                obj.Run();
                PointD nextposition = obj.position; double nextrotation = obj.rotation;
                nextrotation += obj.angularvelocity / fps;
                nextposition.X += obj.speed.X / fps;
                nextposition.Y += obj.speed.Y / fps;
                nextrotation = nextrotation % 360.ToRadians();

                List<Collision> collisions = new List<Collision>();

                foreach (PointD point in obj.hitpoints)
                {
                    //Наступна позиція даної точки
                    PointD col_pos = GetMatrixPosition(point, nextposition, nextrotation);
                    bool? res = GetMatrixState(col_pos);
                    if (res != null)
                    {
                        if (res == true)
                        {
                            double l = 0, r = 1, mid;
                            double res_rotation = 0;
                            PointD res_position = new PointD(0,0);
                            for (int i = 0; i < precision; i++)
                            {
                                mid = (l + r) / 2d;
                                res_rotation = obj.rotation + mid * obj.angularvelocity / fps;
                                res_position = obj.position + mid * obj.speed / fps;
                                if (GetMatrixState(GetMatrixPosition(point, res_position, res_rotation)) == true)
                                {
                                    r = mid;
                                }
                                else
                                {
                                    l = mid;
                                }
                            }
                            PointD speed = (col_pos - GetMatrixPosition(point, res_position, res_rotation)) * fps;
                            double speed_abs = Math.Sqrt(speed.X * speed.X + speed.Y * speed.Y);
                            double surf_angle = GetSurfaceAngle(col_pos.X * pim);
                            double col_angle = Math.Atan2(speed.X, speed.Y) - surf_angle;
                            collisions.Add(new Collision(point, surf_angle, speed_abs * Math.Sin(col_angle), speed_abs * Math.Cos(col_angle), res_rotation));
                            // alpha = 1/2 * pi - col_angle
                            // beta = col_angle
                            // perpendicular_to_suface  = speed * sin(col_angle)
                            // parallel_to_surface = speed * cos(col_angle)
                        }
                    }
                    else
                    {
                        obj.Kill();
                        return;
                    }
                }
                double resmid = 0;
                if (collisions.Count!=0)
                {
                    double initial_r = 1.05;
                    if (collisions.Count >= 2) 
                        initial_r = 1.1 / Math.Pow(Math.Cos((collisions[0].surface_angle - collisions[collisions.Count - 1].surface_angle) / 2), 2);
                    double l = 0, r = initial_r / (double)collisions.Count(), mid;
                    PointD initial_speed = obj.speed;
                    double initial_angular_velociity = obj.angularvelocity;
                    for (int i = 0; i < precision; i++)
                    {
                        obj.speed = initial_speed;
                        obj.angularvelocity = initial_angular_velociity;
                        mid = (l + r) / 2d;
                        resmid = mid;
                        foreach (Collision collision in collisions)
                        {
                            double extrusion_force = Math.Max(0, collision.perpendicular_speed * obj.mass) * mid;
                            double friction_force = Math.Min(Math.Abs(collision.parallel_speed), Math.Max(0, collision.perpendicular_speed) * obj.friction_coefficient) * obj.mass * mid;
                            obj.ApplyForce(collision.pos.Turn(collision.object_rotation), collision.surface_angle - 0.5 * Math.PI, extrusion_force, i == (precision - 1));
                            obj.ApplyForce(collision.pos.Turn(collision.object_rotation), collision.surface_angle, -1 * collision.parallel_speed.GetSign() * friction_force, i == (precision - 1));
                        }
                        bool allout = true;
                        foreach (Collision collision in collisions)
                        {
                            PointD col_pos = GetMatrixPosition(collision.pos, obj.position + obj.speed / fps, obj.rotation + obj.angularvelocity / fps);
                            if (GetMatrixState(col_pos) == true)
                            {
                                allout = false;
                                break;
                            }
                        }
                        if (allout)
                        {
                            r = mid;
                        }
                        else
                        {
                            l = mid;
                        }
                    }
                    obj.position += obj.speed / fps;
                    obj.rotation += obj.angularvelocity / fps;
                }
                else
                {
                    obj.rotation = nextrotation;
                    obj.position = nextposition;
                }
                // check for death points
                // must be done after detection&reaction on colisions
                foreach (PointD point in obj.deathpoints)
                {
                    if (GetMatrixState(GetMatrixPosition(point, obj.position, obj.rotation)).In(null, true))
                    {
                        obj.Kill();
                        return;
                    }
                }
            }
        }*/

        private void RunObjPhysics(MP_Object obj)
        {
            if (obj.alive)
            {
                obj.forces.Clear();
                obj.Run();
                PointD nextposition = obj.position; double nextrotation = obj.rotation;
                nextrotation += obj.angularvelocity / fps;
                nextposition.X += obj.speed.X / fps;
                nextposition.Y += obj.speed.Y / fps;
                nextrotation = nextrotation % 360.ToRadians();

                List<Collision> collisions = new List<Collision>();
                int curp = 0;
                foreach (PointD point in obj.hitpoints)
                {
                    //Наступна позиція даної точки
                    PointD col_pos = GetMatrixPosition(point, nextposition, nextrotation);
                    bool? res = GetMatrixState(col_pos);
                    if (res != null)
                    {
                        if (res == true)
                        {
                            bool? is_inside = GetMatrixState(GetMatrixPosition(point, obj.position, obj.rotation));//якщо позиція вже всередині
                            obj.ever_collided[curp] = true;
                            obj.now_collided[curp] = true;
                            if (is_inside == true)
                            {
                                double surf_angle = GetSurfaceAngle(col_pos.X * pim); // нахил поверхності
                                PointD pixel_position = col_pos * pim; // піксельна позиція(не в метрах, а в пікселях)
                                double depth = GetMatrixHeight(pixel_position.X) - pixel_position.Y; // глубина занурення у поверхню
                                double perp_s = depth * Math.Sin(surf_angle) / pim * fps; // перпендикулярне занурення у поверхню
                                PointD speed = (col_pos - GetMatrixPosition(point, obj.position, obj.rotation)) * fps; // швидкість даної точки[м/с] - вектор
                                double speed_abs = Math.Sqrt(speed.X * speed.X + speed.Y * speed.Y); // скалярна величина швидкості даної точки[м/с]
                                double col_angle = Math.Atan2(speed.X, speed.Y) - surf_angle; // кут зіткнення з поверхнею
                                collisions.Add(new Collision(point, surf_angle, perp_s, speed_abs * Math.Cos(col_angle), nextrotation));
                            }
                            else if (is_inside == false)
                            {
                                double l = 0, r = 1, mid;
                                double res_rotation = 0;
                                PointD res_position = new PointD(0, 0);
                                for (int i = 0; i < precision; i++)
                                {
                                    mid = (l + r) / 2d;
                                    res_rotation = obj.rotation + mid * obj.angularvelocity / fps;
                                    res_position = obj.position + mid * obj.speed / fps;
                                    if (GetMatrixState(GetMatrixPosition(point, res_position, res_rotation)) == true)
                                    {
                                        r = mid;
                                    }
                                    else
                                    {
                                        l = mid;
                                    }
                                }
                                PointD speed = (col_pos - GetMatrixPosition(point, res_position, res_rotation)) * fps;
                                double speed_abs = Math.Sqrt(speed.X * speed.X + speed.Y * speed.Y);
                                double surf_angle = GetSurfaceAngle(col_pos.X * pim);
                                double col_angle = Math.Atan2(speed.X, speed.Y) - surf_angle;
                                collisions.Add(new Collision(point, surf_angle, speed_abs * Math.Sin(col_angle), speed_abs * Math.Cos(col_angle), res_rotation));
                                // alpha = 1/2 * pi - col_angle
                                // beta = col_angle
                                // perpendicular_to_suface  = speed * sin(col_angle)
                                // parallel_to_surface = speed * cos(col_angle)
                            }
                        }
                        else
                        {
                            obj.now_collided[curp] = false;
                        }
                        curp++;
                    }
                    else
                    {
                        obj.Kill();
                        return;
                    }
                }
                double resmid = 0;
                if (collisions.Count != 0)
                {
                    obj.Collide();
                    double initial_r = 1.05;
                    if (collisions.Count >= 2)
                        initial_r = 1.1 / Math.Pow(Math.Cos((collisions[0].surface_angle - collisions[collisions.Count - 1].surface_angle) / 2), 2);
                    double l = 0, r = initial_r / (double)collisions.Count(), mid;
                    PointD initial_speed = obj.speed;
                    double initial_angular_velociity = obj.angularvelocity;
                    for (int i = 0; i < precision; i++)
                    {
                        obj.speed = initial_speed;
                        obj.angularvelocity = initial_angular_velociity;
                        mid = (l + r) / 2d;
                        resmid = mid;
                        foreach (Collision collision in collisions)
                        {
                            double extrusion_force = Math.Max(0, collision.perpendicular_speed * obj.mass) * mid;
                            double friction_force = Math.Min(Math.Abs(collision.parallel_speed), Math.Max(0, collision.perpendicular_speed) * obj.friction_coefficient) * obj.mass * mid;
                            obj.ApplyForce(collision.pos.Turn(collision.object_rotation), collision.surface_angle - 0.5 * Math.PI, extrusion_force, i == (precision - 1));
                            obj.ApplyForce(collision.pos.Turn(collision.object_rotation), collision.surface_angle, -1 * collision.parallel_speed.GetSign() * friction_force, i == (precision - 1));
                        }
                        bool allout = true;
                        foreach (Collision collision in collisions)
                        {
                            PointD col_pos = GetMatrixPosition(collision.pos, obj.position + obj.speed / fps, obj.rotation + obj.angularvelocity / fps);
                            if (GetMatrixState(col_pos) == true)
                            {
                                allout = false;
                                break;
                            }
                        }
                        if (allout)
                        {
                            r = mid;
                        }
                        else
                        {
                            l = mid;
                        }
                    }
                    obj.position += obj.speed / fps;
                    obj.rotation += obj.angularvelocity / fps;
                }
                else
                {
                    obj.rotation = nextrotation;
                    obj.position = nextposition;
                }
                // check for death points
                // must be done after detection&reaction on colisions
                foreach (PointD point in obj.deathpoints)
                {
                    if (GetMatrixState(GetMatrixPosition(point, obj.position, obj.rotation)).In(null, true))
                    {
                        obj.Kill();
                        return;
                    }
                }
            }
        }

        public void Run()
        {
            foreach (MP_Object obj in objects)
            {
                obj.BeforeRun();
                if (obj.active) RunObjPhysics(obj);
                obj.AfterRun();
            }
        }
    }
}
