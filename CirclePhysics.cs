using Microsoft.Xna.Framework;
using MonoHelper;
using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PhysicsLib
{
    public class CircleObject
    {
        public PointD position;
        public double mass;
        public List<Point> myrects = new List<Point>(4);
        public PointD speed;
        public double radius;
        internal bool to_delete = false;

        public CircleObject (PointD _pos, double _mass, double _radius)
        {
            position = _pos;
            mass = _mass;    
            radius = _radius;
        }

        public virtual void Delete()
        {
            to_delete = true;
        }

        public virtual void Run()
        {

        }
    }

    public class CirclePhysics
    {
        public List<CircleObject> circleobjects = new List<CircleObject>();
        public List<int>[,] presens_array;
        public Rectangle map_range;
        int max_radius_of_particle;
        int fps { get { return screen_fps * precision; } }
        public int screen_fps = 0;
        public int precision = 0;

        public CirclePhysics(int _max_radius_of_particle, Rectangle _map_range, int _screen_fps, int _precision)
        {
            max_radius_of_particle = _max_radius_of_particle;
            map_range = _map_range;
            presens_array = new List<int>[_map_range.Width / (max_radius_of_particle * 2) + 4, _map_range.Height / (max_radius_of_particle * 2) + 4];
            for (int i = 0; i < presens_array.GetLength(0); i++)
            {
                for (int j = 0; j < presens_array.GetLength(1); j++)
                {
                    presens_array[i, j] = new List<int>();
                }
            }

            screen_fps = _screen_fps;
            precision = _precision;
        }


        /// <summary>
        /// Returns place of POSITION on PRESENS_ARRAY
        /// </summary>
        /// <param name="position"></param>
        /// <returns></returns>
        public Point Get_rect(PointD position)
        {
            return new Point((int)((position.X - map_range.X) / max_radius_of_particle / 2) + 1, (int)((position.Y - map_range.Y) / max_radius_of_particle / 2) + 1);
        }

        public List<int> Elements_at_rect(Point rect)
        {
            return presens_array[rect.X, rect.Y];
        }

        /// <summary>
        /// Returns index of element located at POSITION. Note that function returns index of only one element.
        /// </summary>
        /// <param name="position"></param>
        /// <returns></returns>
        public int Index_at_position(PointD position)
        {
            if (position.InRect(map_range))
            {
                List<int> indices = Elements_at_rect(Get_rect(position));
                foreach (var index in indices)
                {
                    if ((position - circleobjects[index].position).Length() < circleobjects[index].radius) return index;
                }
            }
            return -1;
        }

        private void Push_index(Point rect_point, int index)
        {
            circleobjects[index].myrects.Add(rect_point);
            presens_array[rect_point.X, rect_point.Y].Add(index);
        }

        List<Point> Get_rects(PointD pos, double radius)
        {
            List<Point> rects = new List<Point>(4);
            Point rect_pos = Get_rect(pos);
            rects.Add(rect_pos);
            int x = 0, y = 0;
            if ((pos.X + radius) / (max_radius_of_particle * 2) != pos.X / (max_radius_of_particle * 2)) x = 1;
            if ((pos.X - radius) / (max_radius_of_particle * 2) != pos.X / (max_radius_of_particle * 2)) x = -1;
            if ((pos.Y + radius) / (max_radius_of_particle * 2) != pos.Y / (max_radius_of_particle * 2)) y = 1;
            if ((pos.Y - radius) / (max_radius_of_particle * 2) != pos.Y / (max_radius_of_particle * 2)) y = -1;
            if (x != 0) rects.Add(new Point(rect_pos.X + x, rect_pos.Y));
            if (y != 0) rects.Add(new Point(rect_pos.X, rect_pos.Y + y));
            if (x != 0 && y != 0) rects.Add(new Point(rect_pos.X + x, rect_pos.Y + y));
            return rects;
        }

        public bool Is_free(PointD pos, double radius)
        {
            var myrects = Get_rects(pos, radius);
            foreach (var rect in myrects)
            {
                if (presens_array[rect.X, rect.Y].Count == 1) continue;
                for (int c = 0; c < presens_array[rect.X, rect.Y].Count; c++)
                {
                    int j = presens_array[rect.X, rect.Y][c];
                    if (j < circleobjects.Count && (pos - circleobjects[j].position).Length() < (circleobjects[j].radius + radius)) return false;
                }
            }
            return true;
        }

        public void MakeLocation()
        {
            for (int i = 0; i < presens_array.GetLength(0); i++)
            {
                for (int j = 0; j < presens_array.GetLength(1); j++)
                {
                    presens_array[i, j].Clear();
                }
            }
            for (int i = 0; i < circleobjects.Count; i++)
            {
                circleobjects[i].myrects.Clear();
                var new_rects = Get_rects(circleobjects[i].position, circleobjects[i].radius);
                foreach (var t_rect in new_rects)
                {
                    Push_index(t_rect, i);
                }
            }
        }

        double FA(double a1, double a2)
        {
            if (a1 < 0) a1 += 2 * Math.PI;
            if (a2 < 0) a2 += 2 * Math.PI;
            double res = Math.Abs(a1 - a2);
            if (res > Math.PI) res = 2 * Math.PI - res;
            return res;
        }

        public virtual void Collide(int i, int j)
        {
            double t = 0;
            PointD posi = circleobjects[i].position, posj = circleobjects[j].position;
            {
                double b = (posi - posj).Length();
                double r = circleobjects[i].radius + circleobjects[j].radius;
                double alpha = FA((circleobjects[j].speed - circleobjects[i].speed).Angle(), (posj - posi).Angle());
                double beta = Math.Asin(b * Math.Sin(alpha) / r);
                double gamma = Math.PI - alpha - beta;
                double c = Math.Sin(gamma) * r / Math.Sin(alpha);
                t = c / (circleobjects[j].speed + circleobjects[i].speed).Length() * fps;
                if (t>=1)
                {
                        
                }

            }
            if (t >= 0 && t <= 1) 
            {
                posi -= t * circleobjects[i].speed / (double)fps;
                posj -= t * circleobjects[j].speed / (double)fps;
                double anglei = (posj - posi).Angle();
                double anglej = (posi - posj).Angle();
                double alpha = FA(circleobjects[i].speed.Angle(), anglei);
                double beta = FA(circleobjects[j].speed.Angle(), anglej);
                double force_i = Math.Cos(alpha) * circleobjects[i].speed.Length();
                double force_j = Math.Cos(beta) * circleobjects[j].speed.Length();
                double mass_coef = circleobjects[i].mass / (circleobjects[i].mass + circleobjects[j].mass);
                circleobjects[j].speed += mass_coef * new PointD((force_j + force_i) * Math.Sin(anglei), (force_j + force_i) * Math.Cos(anglei));
                circleobjects[i].speed += (1 - mass_coef) * new PointD((force_j + force_i) * Math.Sin(anglej), (force_j + force_i) * Math.Cos(anglej));
                posi += t * circleobjects[i].speed / (double)fps;
                posj += t * circleobjects[j].speed / (double)fps;
                circleobjects[i].position = posi; circleobjects[j].position = posj;
            }
            else
            {
                double d = (circleobjects[i].radius + circleobjects[j].radius) - (circleobjects[i].position - circleobjects[j].position).Length();
                if (d > 0)
                {
                    //PointD prev_posi = posi - circleobjects[i].speed / fps, prev_posj = posj - circleobjects[j].speed / fps;
                    double angle = (circleobjects[i].position - circleobjects[j].position).Angle();
                    double mass_coef = circleobjects[i].mass / (circleobjects[i].mass + circleobjects[j].mass);
                    circleobjects[i].position.X += Math.Sin(angle) * d * (1 - mass_coef);
                    circleobjects[i].position.Y += Math.Cos(angle) * d * (1 - mass_coef);
                    circleobjects[j].position.X -= Math.Sin(angle) * d * mass_coef;
                    circleobjects[j].position.Y -= Math.Cos(angle) * d * mass_coef;
                    //circleobjects[i].speed = circleobjects[i].position - prev_posi;
                    //circleobjects[j].speed = circleobjects[j].position - prev_posj;
                }
            }
        }

        public virtual void RunPhysics()
        {
            foreach (var circle in circleobjects)
            {
                circle.Run();
            }
            for (int x = 0; x < precision; x++)
            {
                MakeLocation();
                List<PointD> prevpos = new List<PointD>();
                foreach (CircleObject obj in circleobjects)
                {
                    prevpos.Add(obj.position);
                    obj.position += obj.speed / fps;
                }
                for (int i = 0; i < circleobjects.Count; i++) if (!circleobjects[i].to_delete)
                {
                    foreach (var rect in circleobjects[i].myrects)
                    {
                        if (presens_array[rect.X, rect.Y].Count == 1) continue;
                        for (int c = 0; c < presens_array[rect.X, rect.Y].Count; c++)
                        {
                            int j = presens_array[rect.X, rect.Y][c];
                            if ((i == j) || circleobjects[j].to_delete) continue;
                            if ((circleobjects[i].position - circleobjects[j].position).Length() < (circleobjects[i].radius + circleobjects[j].radius)) Collide(i, j);
                        }
                    }
                }
                for (int i = 0; i < circleobjects.Count; i++)
                {
                    if (circleobjects[i].to_delete)
                    {
                        circleobjects.RemoveAt(i);
                        i--;
                    }
                }
                MakeLocation();
            }
        }
    }
}
