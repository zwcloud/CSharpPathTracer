using System;
using System.Diagnostics;
using System.IO;
using System.Numerics;
using SixLabors.ImageSharp;
using SixLabors.ImageSharp.PixelFormats;

namespace CSharpPathTracer
{
    public class Program
    {
        private struct Color
        {
            public float r, g, b;

            public Color(float r, float g, float b)
            {
                this.r = r;
                this.g = g;
                this.b = b;
            }

            public static Color operator+(Color x, Color y)
            {
                return new Color(x.r + y.r, x.g + y.g, x.b + y.b);
            }

            public static Color operator*(Color x, Color y)
            {
                return new Color(x.r * y.r, x.g * y.g, x.b * y.b);
            }

            public static Color operator *(Color color, float d)
            {
                return new Color(color.r * d, color.g * d, color.b * d);
            }

            public static Color operator /(Color color, float d)
            {
                return new Color(color.r / d, color.g / d, color.b / d);
            }

            public static readonly Color black = new Color(0,0,0);
        }

        private struct Ray
        {
            public Vector3 Origin;
            public Vector3 Direction;
        }

        private class Sphere
        {
            public string Name { get; }
            public Vector3 Center { get; }
            public float Radius { get; }
            public Material Material { get; }

            public Sphere(string name, float radius, Vector3 center, Material material)
            {
                this.Name = name;
                this.Center = center;
                this.Radius = radius;
                this.Material = material;
            }
        }

        private struct Hit
        {
            public Sphere thingHit { get; set; }
            public Vector3 pointWhereObjWasHit { get; set; }
            public Vector3 normalWhereObjWasHit { get; set; }
        }

        private static class Camera
        {
            private static Vector3 Position { get; } = new Vector3(50, 52, 295.6f);
            private static Vector3 Direction { get; } = Vector3.Normalize(new Vector3(0, -0.042612f, -1));

            //ref: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays

            public static Ray generateRay(int x, int y, int w, int h, float fov)
            {
                var ratio = w * 1.0f / h;
                var PixelCameraX = (2 * (x + 0.5f) / w - 1) * ratio * MathF.Tan(fov / 2);
                var PixelCameraY = (1 - 2 * (y + 0.5f) / h) * MathF.Tan(fov / 2);
                var PixelCameraZ = -1;
                var P_cameraSpace = new Vector3(PixelCameraX, PixelCameraY, PixelCameraZ);

                Matrix4x4 M = Matrix4x4.CreateWorld(Position, Direction, Vector3.UnitY);

                var P = Vector3.Transform(P_cameraSpace, M);

                var O = Vector3.Transform(Vector3.Zero, M);

                var OP = P - O;

                var ray = new Ray
                {
                    Origin = O,
                    Direction = Vector3.Normalize(OP)
                };
                return ray;
            }
        }

        private class Material
        {
            public Color Emittance { get; set; }
            public Color Reflectance { get; set; }

        }

        private static Vector3 Project(Vector3 from, Vector3 onto)
        {
            var u = from;
            var v = onto;
            return Vector3.Dot(v, u) / v.LengthSquared() * v;
        }

        /// <summary>
        /// Get the first intersection of a ray and a sphere.
        /// </summary>
        /// <param name="C">center of the sphere</param>
        /// <param name="r">radius of the sphere</param>
        /// <param name="P">origin of the ray</param>
        /// <param name="d">direction vector of the ray</param>
        /// <param name="k"></param>
        /// <param name="I">first intersection point</param>
        /// <param name="n">sphere surface normal at the point</param>
        /// <returns>whether the ray and sphere are intersected</returns>
        private static bool Intersect(Vector3 C, float r, Vector3 P, Vector3 d, out float k, out Vector3 I, out Vector3 n)
        {
            var PC = C - P;
            if (Vector3.Dot(PC, d) < 0)
            {
                k = 0;
                I = Vector3.Zero;
                n = Vector3.Zero;
                return false;
            }
            var PA = Project(PC, d);
            var A = P + PA;
            var CA = A - C;
            var dCA = CA.Length();
            if (dCA > r)
            {
                k = 0;
                I = Vector3.Zero;
                n = Vector3.Zero;
                return false;
            }
            var dIA = MathF.Sqrt(r * r - dCA * dCA);
            var dPA = PA.Length();
            var dPI = dPA - dIA;
            k = dPI;
            I = P + Vector3.Normalize(d) * dPI;
            var CI = I - C;
            n = CI / r;//normalized
            return true;
        }

        private static Random Random = new Random();

        private static Vector3 RandomUnitVectorOnUnitSphere()
        {
            //ref: http://mathworld.wolfram.com/SpherePointPicking.html
            var x0 = Random.NextDouble()*2 - 1;
            var x1 = Random.NextDouble()*2 - 1;
            var x2 = Random.NextDouble()*2 - 1;
            var x3 = Random.NextDouble()*2 - 1;
            var divider = x0 * x0 + x1 * x1 + x2 * x2 + x3 * x3;
            var pX = (float)(2 * (x1 * x3 + x0 * x2) / divider);
            var pY = (float)(2 * (x2 * x3 - x0 * x1) / divider);
            var pZ = (float)((x0 * x0 + x3 * x3 - x1 * x1 - x2 * x2) / divider);
            return new Vector3(pX, pY, pZ);
        }

        private static Vector3 RandomUnitVectorOnNorthernHemisphere()
        {
            var p = RandomUnitVectorOnUnitSphere();
            Debug.Assert(Math.Abs(p.LengthSquared()-1) <= 0.001);
            //move the point onto northern hemisphere surface
            p.Y = Math.Abs(p.Y);

            return p;
        }

        private static float AngleBetween(Vector3 a, Vector3 b)
        {
            return (float)Math.Acos(Vector3.Dot(a, b) / (a.Length() * b.Length()));
        }

        private static Vector3 RotateUnitVector(Vector3 p, Vector3 a, Vector3 b)
        {
            a = Vector3.Normalize(a);
            b = Vector3.Normalize(b);
            var axis = Vector3.Normalize(Vector3.Cross(b, a));
            var angle = AngleBetween(a, b);
            var quaternion = Quaternion.CreateFromAxisAngle(axis, angle);
            return Vector3.Transform(p, quaternion);
        }

        private static Vector3 RandomUnitVectorInHemisphereOf(Vector3 dir)
        {
            var p = RandomUnitVectorOnNorthernHemisphere();
            //now p is distributed around the north unit vector: Vector3.UnitY
            p = RotateUnitVector(p, Vector3.Normalize(dir), Vector3.UnitY);//rotate the vector to make it surround dir
            //now p is distributed around dir
            Debug.Assert(Math.Abs(p.LengthSquared() - 1) < 0.001);
            return p;
        }

        private static readonly Sphere[] Spheres =
        {
            new Sphere("Sphere0", 500, new Vector3(7,327,-875), new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.25f,0.25f)}),//Wall0
            new Sphere("Sphere1", 500, new Vector3(-769,55,-94),  new Material{Emittance = Color.black, Reflectance = new Color(0.25f,0.25f,0.75f)}),//Wall1
            new Sphere("Sphere2", 500, new Vector3(720,16, -58),  new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.75f,0.75f)}),//Wall2
            new Sphere("Ball",  50,  new Vector3(15,16.5f,-91), new Material{Emittance = Color.black, Reflectance = new Color(0.25f,0.25f,0.75f)}),//Ball
            new Sphere("Light", 50,  new Vector3(-148,52,-60),  new Material{Emittance = new Color(12,12,12), Reflectance = Color.black})       //Light
        };

        private const int MaxDepth = 10;
        private static Color TracePath(Ray ray, int depth)
        {
            if (depth > MaxDepth)
            {
                return Color.black;// Bounced enough times.
            }

            if (!FindNearestObject(ray, out var hit))
            {
                return Color.black;// Nothing was hit.
            }

            Material material = hit.thingHit.Material;
            Color emittance = material.Emittance;

            // Pick a random direction from here and keep going.
            Ray newRay;
            newRay.Origin = hit.pointWhereObjWasHit;

            // This is NOT a cosine-weighted distribution!
            newRay.Direction = RandomUnitVectorInHemisphereOf(hit.normalWhereObjWasHit);

            // Probability of the newRay
            const float p = 1/(2*MathF.PI);

            // Compute the BRDF for this ray (assuming Lambertian reflection)
            float cos_theta = Vector3.Dot(newRay.Direction, hit.normalWhereObjWasHit);
            Color BRDF = material.Reflectance / MathF.PI;

            // Recursively trace reflected light sources.
            Color incoming = TracePath(newRay, depth + 1);

            // Apply the Rendering Equation here.
            return emittance + (BRDF * incoming * cos_theta / p);
        }

        private static bool FindNearestObject(Ray ray, out Hit hit)
        {
            hit = new Hit();
            float minK = float.MaxValue;
            foreach (var sphere in Spheres)
            {
                float k;
                bool isHit = Intersect(sphere.Center, sphere.Radius, ray.Origin, ray.Direction, out k,
                    out var I, out var n);
                if (isHit)
                {
                    if (k < minK)
                    {
                        minK = k;
                        hit.thingHit = sphere;
                        hit.pointWhereObjWasHit = I;
                        hit.normalWhereObjWasHit = n;
                    }
                }
            }
            if (minK != float.MaxValue)
            {
                return true;
            }
            return false;
        }

        private static void Render(Color[] image, int width, int height, int numSamples)
        {
            const float fov = MathF.PI / 3;
            for (var y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    var index = y * width + x;
                    var pixel = image[index];
                    Ray r = Camera.generateRay(x, y, width, height, fov);
                    for (var i = 0; i < numSamples; i++)
                    {
                        pixel += TracePath(r, 0);
                    }
                    pixel /= numSamples;
                    image[index] = pixel;
                }
            }
        }

        public static void Main(string[] args)
        {
            const int w = 400, h = 300;
            var finalImage = new Color[w * h];
            var numSamples = 10;
            Render(finalImage, w, h, numSamples);

            //save image
            var bytes = new byte[w * h * 3];
            for (var i = 0; i < finalImage.Length; i++)
            {
                var color = finalImage[i];
                byte r = (byte) (color.r * 255);
                byte g = (byte) (color.g * 255);
                byte b = (byte) (color.b * 255);
                bytes[3*i] = r;
                bytes[3*i+1] = g;
                bytes[3*i+2] = b;
            }
            using (var image = Image.LoadPixelData<Rgb24>(bytes, w, h))
            using (var stream = new FileStream("D:\\1.png", FileMode.OpenOrCreate))
            {
                image[0,0] = new Rgb24(255, 0, 0);
                image[image.Width-1,image.Height-1] = new Rgb24(0, 255, 0);
                image.SaveAsPng(stream);
            }
        }

    }
}
