# CSharpPathTracer

A naive [path tracing](https://en.wikipedia.org/wiki/Path_tracing) implementation written in C#.

## Build & Run

Open the solution in Visual Studio 2017, then build and run.

Or build and run with .NET Core:

```plain
cd CSharpPathTracer
dotnet run
```

Then the rendered image is saved to
* Windows: `C:\Users\<your_user_name>\xxx.png`
* Linux/macOS: `~/xxx.png`.

## Result

Environment: Win10, i7-8700

Sample: 4096, Depth: 10, Time Used: 641.88s  
![](Images/4096_641.88.png)

## Article

* (In Chinese) [Path Tracing](https://zwcloud.net/#blog/98)
* (In Chinese) [BRDF](https://zwcloud.net/#blog/99)
* (In Chinese) [Material: Ideal Diffuse and Ideal Specular](https://zwcloud.net/#blog/100)
* (In Chinese) [Refraction: Glass Material](https://zwcloud.net/#blog/101)
