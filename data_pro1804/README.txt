BEM input files used in evaluations of the PRO-2018-4 paper.

input_50ms.bin:
- Sphere: a sphere composed of 50,000,000 particles.
  The number of particles = 50,000,000

input_100ts.txt_cb_1.5_10_10_10.bin:
- SphereCube: $10\times 10\times 10$ spheres placed cubically. Each sphere is composed of 101,250 particles.
  The number of particles = 101,250,000

input_100ts.txt_pb_1.5_14.bin:
- SpherePyramid: $1^2+2+2+\ldots+14^2=1015$ spheres placed pyramidally. Each sphere is composed of 101,250 particles.
  The number of particles = 102,768,750

input_human_1x1.txt_cb_0.3_50_100_1.bin:
- Humans: $50\times 100$ pairs of human-shaped objects. Each pair of objects is composed of 19,664 particles.
* The picture (input_human_1x1.txt_cb_0.3_5_10_1.png) shows only $5\times 10$ pairs of objects.

The faces of sphere and human-shaped objects are composed of triangles. The gravity centers of the triangles
are treated as particles in our evaluations.
