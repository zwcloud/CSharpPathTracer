using System;
using System.Diagnostics;
using System.IO;
using System.Numerics;
using System.Threading;
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

        enum MaterialType
        {
            Diffuse,
            Mirror,
            Glass,
        }

        private class Material
        {
            public Color Emittance { get; set; }
            public Color Reflectance { get; set; }
            public MaterialType Type { get; set; }
        }

        private static Vector3 Project(Vector3 from, Vector3 onto)
        {
            var u = from;
            var v = onto;
            return Vector3.Dot(v, u) / v.LengthSquared() * v;
        }

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
            var dPC = PC.Length();
            var dPI = (dPC-r>0.001) ? (dPA - dIA) : (dPA + dIA);
            k = dPI;
            I = P + Vector3.Normalize(d) * dPI;
            var CI = I - C;
            n = CI / r;//normalized
            return true;
        }

        //ref: https://codeblog.jonskeet.uk/2009/11/04/revisiting-randomness/
        public static class ThreadLocalRandom
        {
            private static readonly Random globalRandom = new Random();
            private static readonly object globalLock = new object();

            private static readonly ThreadLocal<Random> threadRandom = new ThreadLocal<Random>(NewRandom);

            public static Random NewRandom()
            {
                lock (globalLock)
                {
                    return new Random(globalRandom.Next());
                }
            }

            public static Random Instance { get { return threadRandom.Value; } }

            public static int Next()
            {
                return Instance.Next();
            }
        }

        public static Vector3 RandomUnitVectorOnUnitSphere()
        {
            //ref: http://mathworld.wolfram.com/SpherePointPicking.html
            var x0 = ThreadLocalRandom.Instance.NextDouble()*2 - 1;
            var x1 = ThreadLocalRandom.Instance.NextDouble()*2 - 1;
            var x2 = ThreadLocalRandom.Instance.NextDouble()*2 - 1;
            var x3 = ThreadLocalRandom.Instance.NextDouble()*2 - 1;
            var divider = x0 * x0 + x1 * x1 + x2 * x2 + x3 * x3;
            var pX = (float)(2 * (x1 * x3 + x0 * x2) / divider);
            var pY = (float)(2 * (x2 * x3 - x0 * x1) / divider);
            var pZ = (float)((x0 * x0 + x3 * x3 - x1 * x1 - x2 * x2) / divider);
            return new Vector3(pX, pY, pZ);
        }

        public static Vector3 RandomUnitVectorOnNorthernHemisphere()
        {
            var p = RandomUnitVectorOnUnitSphere();
            Debug.Assert(Math.Abs(p.LengthSquared()-1) <= 0.001);
            //move the point onto northern hemisphere surface
            p.Y = Math.Abs(p.Y);

            return p;
        }

        public static float AngleBetween(Vector3 a, Vector3 b)
        {
            return (float)Math.Acos(Vector3.Dot(a, b) / (a.Length() * b.Length()));
        }

        public static Vector3 RotateUnitVector(Vector3 p, Vector3 a, Vector3 b)
        {
            a = Vector3.Normalize(a);
            b = Vector3.Normalize(b);
            var axis = Vector3.Normalize(Vector3.Cross(a, b));
            var angle = AngleBetween(a, b);
            var quaternion = Quaternion.CreateFromAxisAngle(axis, angle);
            return Vector3.Transform(p, quaternion);
        }

        public static Vector3 RandomUnitVectorInHemisphereOf(Vector3 dir)
        {
            var p = RandomUnitVectorOnNorthernHemisphere();
            //now p is distributed around the north unit vector: Vector3.UnitY
            p = RotateUnitVector(p, Vector3.UnitY, Vector3.Normalize(dir));//rotate the vector to make it surround dir
            //now p is distributed around dir
            return p;
        }

        private static readonly Sphere[] Spheres =
        {
            new Sphere("Sphere0", 500, new Vector3(7,327,-875), new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.25f,0.25f), Type = MaterialType.Diffuse}),
            new Sphere("Sphere1", 500, new Vector3(-769,55,-94),  new Material{Emittance = Color.black, Reflectance = new Color(0.25f,0.25f,0.75f), Type = MaterialType.Diffuse}),
            new Sphere("Sphere2", 500, new Vector3(720,16, -58),  new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.75f,0.75f), Type = MaterialType.Diffuse}),
            new Sphere("MirrorBall",  50,  new Vector3(15,16.5f,-91), new Material{Emittance = Color.black, Reflectance = new Color(1,1,1), Type = MaterialType.Mirror}),
            new Sphere("GlassBall", 100,  new Vector3(-80,-197,-317),  new Material{Emittance = Color.black, Reflectance = new Color(1,1,1), Type = MaterialType.Glass}),
            new Sphere("LightBall", 50,  new Vector3(-148,52,-60),  new Material{Emittance = new Color(1,1,1)*12, Reflectance = Color.black, Type = MaterialType.Diffuse}),

            new Sphere("Front", 5000, new Vector3(368,16.5f,-5403), new Material{Emittance = Color.black, Reflectance = new Color(0,0.75f,0), Type = MaterialType.Diffuse}),
            new Sphere("Back", 5000, new Vector3(15,16.5f,5073), new Material{Emittance = Color.black, Reflectance = Color.black, Type = MaterialType.Diffuse}),
            new Sphere("Top", 5000, new Vector3(15,5326,91), new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.75f,0.75f), Type = MaterialType.Diffuse}),
            new Sphere("Bottom", 5000, new Vector3(15,-5261,91), new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.75f,0.75f), Type = MaterialType.Diffuse}),
            new Sphere("Left", 5000, new Vector3(-4789,16.5f,-2859), new Material{Emittance = Color.black, Reflectance = new Color(0.75f,0.25f,0.25f), Type = MaterialType.Diffuse}),
            new Sphere("Right", 5000, new Vector3(5487,16.5f,-1256), new Material{Emittance = Color.black, Reflectance = new Color(0.25f,0.25f,0.75f), Type = MaterialType.Diffuse}),
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

            Debug.Assert(
                Math.Abs(hit.pointWhereObjWasHit.X - ray.Origin.X) > 0.0001f
                || Math.Abs(hit.pointWhereObjWasHit.Y - ray.Origin.Y) > 0.0001f
                || Math.Abs(hit.pointWhereObjWasHit.Z - ray.Origin.Z) > 0.0001f
                );

            Material material = hit.thingHit.Material;

            if (material.Type == MaterialType.Diffuse)
            {
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
            else if (material.Type == MaterialType.Mirror)
            {
                Color emittance = material.Emittance;

                // 计算反射光方向
                // 从碰撞点处依照反射定律，发射新的射线
                Ray newRay;
                newRay.Origin = hit.pointWhereObjWasHit;
                // r = i−2(i*n)n
                newRay.Direction = Vector3.Normalize(ray.Direction - 2 * Vector3.Dot(ray.Direction, hit.normalWhereObjWasHit) * hit.normalWhereObjWasHit);

                // 计算反射光颜色
                // 即相同的颜色
                Color BRDF = material.Reflectance;

                // 递归地追踪反射光的来源
                Color incoming = TracePath(newRay, depth + 1);

                // 在此应用渲染方程，得到出射光的最终颜色（强度）
                return emittance + BRDF * incoming;
            }
            else if (material.Type == MaterialType.Glass)
            {
                Color emittance = material.Emittance;

                // 已知
                // 碰撞点位置`ray.pointWhereObjWasHit`
                // 碰撞点处的表面法线`ray.normalWhereObjWasHit`
                // 空气的折射率1，玻璃的折射率1.5
                const float n_air = 1;
                const float n_glass = 1.5f;

                // 判断是从空气入射玻璃还是从玻璃入射空气
                bool intoGlass = Vector3.Dot(ray.Direction, hit.normalWhereObjWasHit) < 0;
                // 计算时的法线方向要根据入射方是否是玻璃球判断是否需要取反
                Vector3 n = intoGlass ? hit.normalWhereObjWasHit : -hit.normalWhereObjWasHit;
                Vector3 d = ray.Direction;

                // 计算反射光方向
                // 从碰撞点处依照反射定律，发射新的反射射线
                Ray reflectionRay;
                reflectionRay.Origin = hit.pointWhereObjWasHit;
                // r = i−2(i*n)n
                reflectionRay.Direction = Vector3.Normalize(d - 2*Vector3.Dot(d, n)*n);

                float n1 = intoGlass ? n_air : n_glass;//射入方的折射率
                float n2 = intoGlass ? n_glass : n_air;//射出方的折射率

                if(!intoGlass)
                {
                    // 判断是否发生全内反射
                    float a = Vector3.Dot(n, d);
	                float b = n2/n1;
                    bool totalInternalRefelctionHappened = 1 > a*a + b*b;
	                if(totalInternalRefelctionHappened)
	                {
		                // 递归地追踪反射光的来源
		                Color incoming = TracePath(reflectionRay, depth + 1);
		                return emittance + material.Reflectance * incoming;
	                }
                }

                // 计算反射光、折射光强度//√
                float R0 = ((n1-n2)/(n1+n2))*((n1-n2)/(n1+n2));//垂直入射时的反射光比例
                float cosineTheta1 = Vector3.Dot(n, -d);
                float R = R0+(1-R0)*MathF.Pow(1-cosineTheta1, 5);//反射光比例，依据菲涅尔方程的估计公式
                float T = 1-R;//折射光比例

                if (!(0<=R && R<=1))
                {
                    Debugger.Break();
                }

                // 计算折射光方向//√
                Ray refractionRay;
                refractionRay.Origin = hit.pointWhereObjWasHit;
                float nd = Vector3.Dot(n, d);
                refractionRay.Direction = Vector3.Normalize(n1*(d-nd*n)/n2 - n*MathF.Sqrt(1-n1*n1*(1-nd*nd)/(n2*n2)));

                // 递归地追踪反射光、折射光的来源
                Color reflection = TracePath(reflectionRay, depth + 1) * R;
                Color refraction = TracePath(refractionRay, depth + 1) * T;

                // 应用渲染方程，得到光的最终颜色（强度）
                return emittance + material.Reflectance * (reflection + refraction);
            }

            throw new ArgumentOutOfRangeException();
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

            System.Threading.Tasks.Parallel.For(0, width*height, index =>
            {
                var x = index % width;
                var y = index / width;
                Debug.Assert(index == y * width + x);

                var pixel = image[index];
                Ray r = Camera.generateRay(x, y, width, height, fov);
                for (var i = 0; i < numSamples; i++)
                {
                    pixel += TracePath(r, 0);
                }
                pixel /= numSamples;
                image[index] = pixel;
            });
        }

        public static void Main(string[] args)
        {
            const int w = 400, h = 300;
            var finalImage = new Color[w * h];
            var numSamples = 4096;

            var watch = Stopwatch.StartNew();
            Render(finalImage, w, h, numSamples);
            var timeUsed = watch.Elapsed.TotalSeconds;

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
            using (var stream = new FileStream($"{Environment.GetFolderPath(Environment.SpecialFolder.UserProfile)}/{numSamples}_{timeUsed:F2}.png", FileMode.OpenOrCreate))
            {
                image[0,0] = new Rgb24(255, 0, 0);
                image[image.Width-1,image.Height-1] = new Rgb24(0, 255, 0);
                image.SaveAsPng(stream);
            }
        }

    }
}
