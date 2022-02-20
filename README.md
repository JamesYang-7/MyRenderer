# MyRenderer

[![license: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)


## Table of Contents

- [MyRenderer](#myrenderer)
  - [Table of Contents](#table-of-contents)
  - [Background](#background)
  - [My Extensions](#my-extensions)
  - [Build](#build)
      - [VS 2019 on Windows](#vs-2019-on-windows)
  - [Usage](#usage)
  - [Examples](#examples)
  - [License](#license)

## Background

This is my exercise program built for learning ray tracing. The code is mainly based on [Ray Tracing in One Weekend Book Series](https://github.com/RayTracing/raytracing.github.io). 

## My Extensions

- Optimized file structure, which makes it easier to add customized functions.
- Photon mapping rendering method.

## Build

```
$ git clone https://github.com/JamesYang-7/MyRenderer.git
```

#### VS 2019 on Windows
1. Open forder `Myrenderer` with VS 2019. 
2. Create a project from existing code.
   [How to create a project from existing code](https://docs.microsoft.com/en-us/cpp/build/how-to-create-a-cpp-project-from-existing-code?view=msvc-170)
3. Open file `project-name.sln` with VS 2019.
4. On the `Build` menu, click `Build Solution`.

## Usage

This is a command line program. In VS 2019, you can click `Debug > Start Without Debugging` on the menu. There are some fixed choice for scenes. In general, you can choose `path tracing` rendering method and `default` settings. The output file format is `.ppm`, you can open it with Photoshop or convert it to general image file formats using some tools.

## Examples

![Cornell box](https://github.com/JamesYang-7/helloworld7/blob/master/images/MyRenderer/cornell_box_hybrid_rate40_50spp.png?raw=true)
![Cornell box with spheres - ray tracing](https://github.com/JamesYang-7/helloworld7/blob/master/images/MyRenderer/cornell_spheres_pt_256spp_32maxdep.png?raw=true)
![Cornell box with spheres - photon mapping](https://github.com/JamesYang-7/helloworld7/blob/master/images/MyRenderer/total_scene.png?raw=true)
## License

[MIT © James Yang.](./LICENSE)

